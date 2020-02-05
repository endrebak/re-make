(ns re-make.core
  "Official Sente reference example: server"
  {:author "Peter Taoussanis (@ptaoussanis)"}

  (:require
   [clojure.java.io :as io]
   [re-make.layout :as layout]
   [re-make.middleware :as middleware]
   [clojure.string     :as str]
   [ring.middleware.defaults]
   [ring.middleware.reload :refer [wrap-reload]]
   [ring.middleware.anti-forgery :as anti-forgery]
   [ring.middleware.webjars :refer [wrap-webjars]]
   [ring.util.http-response :as response]
   [compojure.core     :as comp :refer (defroutes GET POST)]
   [compojure.route    :as route]
   [hiccup.core        :as hiccup]
   [clojure.core.async :as async  :refer (<! <!! >! >!! put! chan go go-loop)]
   [taoensso.encore    :as encore :refer (have have?)]
   [taoensso.timbre    :as timbre :refer (tracef debugf infof warnf errorf)]
   [taoensso.sente     :as sente]

   ;;; TODO Choose (uncomment) a supported web server + adapter -------------
   ;; [org.httpkit.server :as http-kit]
   ;; [taoensso.sente.server-adapters.http-kit :refer (get-sch-adapter)]
   ;;
   [immutant.web :as immutant]
   [taoensso.sente.server-adapters.immutant :refer (get-sch-adapter)]
   ;;
   ;; [nginx.clojure.embed :as nginx-clojure]
   ;; [taoensso.sente.server-adapters.nginx-clojure :refer (get-sch-adapter)]
   ;;
   ;; [aleph.http :as aleph]
   ;; [taoensso.sente.server-adapters.aleph :refer (get-sch-adapter)]
   ;; -----------------------------------------------------------------------

   ;; Optional, for Transit encoding:
   [taoensso.sente.packers.transit :as sente-transit]))

;; (timbre/set-level! :trace) ; Uncomment for more logging
(reset! sente/debug-mode?_ true) ; Uncomment for extra debug info

;;;; Define our Sente channel socket (chsk) server

(let [;; Serializtion format, must use same val for client + server:
      packer :edn ; Default packer, a good choice in most cases
      ;; (sente-transit/get-transit-packer) ; Needs Transit dep

      chsk-server
      (sente/make-channel-socket-server!
       (get-sch-adapter) {:packer packer})

      {:keys [ch-recv send-fn connected-uids
              ajax-post-fn ajax-get-or-ws-handshake-fn]}
      chsk-server]

  (def ring-ajax-post                ajax-post-fn)
  (def ring-ajax-get-or-ws-handshake ajax-get-or-ws-handshake-fn)
  (def ch-chsk                       ch-recv) ; ChannelSocket's receive channel
  (def chsk-send!                    send-fn) ; ChannelSocket's send API fn
  (def connected-uids                connected-uids) ; Watchable, read-only atom
  )

;; We can watch this atom for changes if we like
(add-watch connected-uids :connected-uids
  (fn [_ _ old new]
    (when (not= old new)
      (infof "Connected uids change: %s" new))))

;;;; Ring handlers


(defn login-handler
  "Here's where you'll add your server-side login/auth procedure (Friend, etc.).
  In our simplified example we'll just always successfully authenticate the user
  with whatever user-id they provided in the auth request."
  [ring-req]
  (let [{:keys [session params]} ring-req
        {:keys [user-id]} params]
    (debugf "Login request: %s" params)
    {:status 200 :session (assoc session :uid user-id)}))

(defn home-page [request]
  (layout/render request "home.html"))

(defroutes ring-routes
  (GET  "/"      ring-req (home-page ring-req))
  (GET  "/docs"  ring-req (fn [_]
                            (-> (response/ok (-> "docs/docs.md" io/resource slurp))
                                (response/header "Content-Type" "text/plain; charset=utf-8"))))
  (GET  "/chsk"  ring-req (ring-ajax-get-or-ws-handshake ring-req))
  (POST "/chsk"  ring-req (ring-ajax-post                ring-req))
  (POST "/login" ring-req (login-handler                 ring-req))
  (route/resources "/") ; Static files, notably public/main.js (our cljs target)
  (route/not-found "<h1>Page not found</h1>"))

(def main-ring-handler
  "**NB**: Sente requires the Ring `wrap-params` + `wrap-keyword-params`
  middleware to work. These are included with
  `ring.middleware.defaults/wrap-defaults` - but you'll need to ensure
  that they're included yourself if you're not using `wrap-defaults`.
  You're also STRONGLY recommended to use `ring.middleware.anti-forgery`
  or something similar."
  (ring.middleware.defaults/wrap-defaults
   (-> ring-routes wrap-webjars wrap-reload) ring.middleware.defaults/site-defaults))

;;;; Some server>user async push examples

(defn test-fast-server>user-pushes
  "Quickly pushes 100 events to all connected users. Note that this'll be
  fast+reliable even over Ajax!"
  []
  (doseq [uid (:any @connected-uids)]
    (doseq [i (range 100)]
      (chsk-send! uid [:fast-push/is-fast (str "hello " i "!!")]))))

(comment (test-fast-server>user-pushes))

(defonce broadcast-enabled?_ (atom true))

(defn start-example-broadcaster!
  "As an example of server>user async pushes, setup a loop to broadcast an
  event to all connected users every 10 seconds"
  []
  (let [broadcast!
        (fn [i]
          (let [uids (:any @connected-uids)]
            (debugf "Broadcasting server>user: %s uids" (count uids))
            (doseq [uid uids]
              (chsk-send! uid
                [:some/broadcast
                 {:what-is-this "An async broadcast pushed from server"
                  :how-often "Every 10 seconds"
                  :to-whom uid
                  :i i}]))))]

    (go-loop [i 0]
      (<! (async/timeout 10000))
      (when @broadcast-enabled?_ (broadcast! i))
      (recur (inc i)))))

;;;; Sente event handlers

(defmulti -event-msg-handler
  "Multimethod to handle Sente `event-msg`s"
  :id ; Dispatch on event-id
  )

(defn event-msg-handler
  "Wraps `-event-msg-handler` with logging, error catching, etc."
  [{:as ev-msg :keys [id ?data event]}]
  (-event-msg-handler ev-msg) ; Handle event-msgs on a single thread
  ;; (future (-event-msg-handler ev-msg)) ; Handle event-msgs on a thread pool
  )

(defmethod -event-msg-handler
  :default ; Default/fallback case (no other matching handler)
  [{:as ev-msg :keys [event id ?data ring-req ?reply-fn send-fn]}]
  (let [session (:session ring-req)
        uid     (:uid     session)]
    (debugf "Unhandled event: %s" event)
    (when ?reply-fn
      (?reply-fn {:umatched-event-as-echoed-from-server event}))))

(defmethod -event-msg-handler :example/test-rapid-push
  [ev-msg] (test-fast-server>user-pushes))

(defmethod -event-msg-handler :example/toggle-broadcast
  [{:as ev-msg :keys [?reply-fn]}]
  (let [loop-enabled? (swap! broadcast-enabled?_ not)]
    (?reply-fn loop-enabled?)))

;; TODO Add your (defmethod -event-msg-handler <event-id> [ev-msg] <body>)s here...

;;;; Sente event router (our `event-msg-handler` loop)

(defonce router_ (atom nil))
(defn  stop-router! [] (when-let [stop-fn @router_] (stop-fn)))
(defn start-router! []
  (stop-router!)
  (reset! router_
    (sente/start-server-chsk-router!
      ch-chsk event-msg-handler)))

;;;; Init stuff

(defonce    web-server_ (atom nil)) ; (fn stop [])
(defn  stop-web-server! [] (when-let [stop-fn @web-server_] (stop-fn)))
(defn start-web-server! [& [port]]
  (stop-web-server!)
  (let [port (or port 0) ; 0 => Choose any available port
        ring-handler (var main-ring-handler)

        [port stop-fn]
        ;;; TODO Choose (uncomment) a supported web server ------------------
        ;; (let [stop-fn (http-kit/run-server ring-handler {:port port})]
        ;;   [(:local-port (meta stop-fn)) (fn [] (stop-fn :timeout 100))])
        ;;
        (let [server (immutant/run ring-handler :port port)]
          [(:port server) (fn [] (immutant/stop server))])
        ;;
        ;; (let [port (nginx-clojure/run-server ring-handler {:port port})]
        ;;   [port (fn [] (nginx-clojure/stop-server))])
        ;;
        ;; (let [server (aleph/start-server ring-handler {:port port})
        ;;       p (promise)]
        ;;   (future @p) ; Workaround for Ref. https://goo.gl/kLvced
        ;;   ;; (aleph.netty/wait-for-close server)
        ;;   [(aleph.netty/port server)
        ;;    (fn [] (.close ^java.io.Closeable server) (deliver p nil))])
        ;; ------------------------------------------------------------------

        uri (format "http://localhost:%s/" port)]

    (infof "Web server is running at `%s`" uri)
    (try
      (.browse (java.awt.Desktop/getDesktop) (java.net.URI. uri))
      (catch java.awt.HeadlessException _))

    (reset! web-server_ stop-fn)))

(defn stop!  []  (stop-router!)  (stop-web-server!))
(defn start! [] (start-router!) (start-web-server!) (start-example-broadcaster!))

(defn -main "For `lein run`, etc." [] (start!))

(comment
  (start!)
  (test-fast-server>user-pushes))


;; (ns re-make.core
;;   (:require
;;     [re-make.handler :as handler]
;;     [re-make.nrepl :as nrepl]
;;     [luminus.http-server :as http]
;;     [re-make.config :refer [env]]
;;     [clojure.tools.cli :refer [parse-opts]]
;;     [clojure.tools.logging :as log]
;;     [mount.core :as mount])
;;   (:gen-class))

;; ;; log uncaught exceptions in threads
;; (Thread/setDefaultUncaughtExceptionHandler
;;   (reify Thread$UncaughtExceptionHandler
;;     (uncaughtException [_ thread ex]
;;       (log/error {:what :uncaught-exception
;;                   :exception ex
;;                   :where (str "Uncaught exception on" (.getName thread))}))))

;; (def cli-options
;;   [["-p" "--port PORT" "Port number"
;;     :parse-fn #(Integer/parseInt %)]])

;; (mount/defstate ^{:on-reload :noop} http-server
;;   :start
;;   (http/start
;;     (-> env
;;         (assoc  :handler (handler/app))
;;         (update :io-threads #(or % (* 2 (.availableProcessors (Runtime/getRuntime)))))
;;         (update :port #(or (-> env :options :port) %))))
;;   :stop
;;   (http/stop http-server))

;; (mount/defstate ^{:on-reload :noop} repl-server
;;   :start
;;   (when (env :nrepl-port)
;;     (nrepl/start {:bind (env :nrepl-bind)
;;                   :port (env :nrepl-port)}))
;;   :stop
;;   (when repl-server
;;     (nrepl/stop repl-server)))


;; (defn stop-app []
;;   (doseq [component (:stopped (mount/stop))]
;;     (log/info component "stopped"))
;;   (shutdown-agents))

;; (defn start-app [args]
;;   (doseq [component (-> args
;;                         (parse-opts cli-options)
;;                         mount/start-with-args
;;                         :started)]
;;     (log/info component "started"))
;;   (.addShutdownHook (Runtime/getRuntime) (Thread. stop-app)))

;; (defn -main [& args]
;;   (start-app args))

;; (defn landing-pg-handler [ring-req]
;;   (hiccup/html
;;     [:h1 "Sente reference example"]
;;     (let [csrf-token
;;           ;; (:anti-forgery-token ring-req) ; Also an option
;;           (force anti-forgery/*anti-forgery-token*)]

;;       [:div#sente-csrf-token {:data-csrf-token csrf-token}])
;;     [:p "An Ajax/WebSocket" [:strong " (random choice!)"] " has been configured for this example"]
;;     [:hr]
;;     [:p [:strong "Step 1: "] " try hitting the buttons:"]
;;     [:p
;;      [:button#btn1 {:type "button"} "chsk-send! (w/o reply)"]]
;;      ;; [:button#btn2 {:type "button"} "chsk-send! (with reply)"]]
;;     [:p
;;      [:button#btn3 {:type "button"} "Test rapid server>user async pushes"]
;;      [:button#btn4 {:type "button"} "Toggle server>user async broadcast push loop"]]
;;     [:p
;;      [:button#btn5 {:type "button"} "Disconnect"]
;;      [:button#btn6 {:type "button"} "Reconnect"]]
;;     ;;
;;     [:p [:strong "Step 2: "] " observe std-out (for server output) and below (for client output):"]
;;     [:textarea#output {:style "width: 100%; height: 200px;"}]
;;     ;;
;;     [:hr]
;;     [:h2 "Step 3: try login with a user-id"]
;;     [:p  "The server can use this id to send events to *you* specifically."]
;;     [:p
;;      [:input#input-login {:type :text :placeholder "User-id"}]
;;      [:button#btn-login {:type "button"} "Secure login!"]]
;;     ;;
;;     [:hr]
;;     [:h2 "Step 4: want to re-randomize Ajax/WebSocket connection type?"]
;;     [:p "Hit your browser's reload/refresh button"]
;;     [:div#content]
;;     [:script "js/app.js"] ; Include our cljs target
;;     ))
