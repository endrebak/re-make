(ns re-make.core
  "Official Sente reference example: server"
  {:author "Peter Taoussanis (@ptaoussanis)"}

  (:require
   [clojure.java.io :as io]
   [re-make.watch :as watch]
   [nrepl.server :refer [start-server stop-server]]
   [re-make.config :refer [env]]
   [cprop.source :as source]
   [mount.core :as mount]
   [re-make.layout :as layout]
   [re-make.middleware :as middleware]
   [clojure.string     :as str]
   [ring.middleware.defaults]
   [ring.middleware.reload :refer [wrap-reload]]
   [ring.middleware.anti-forgery :as anti-forgery]
   [ring.middleware.webjars :refer [wrap-webjars]]
   [ring.middleware.format :refer [wrap-restful-format]]
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

(defonce nrepl-server (start-server :bind "127.0.0.1" :port 7000))


;; (println env)
;; (println (source/from-system-props))
;; (println (:nrepl-port (source/from-system-props)))

;; (println (source/from-env))
;; (println (:nrepl-port (source/from-env)))

;;;; Define our Sente channel socket (chsk) server

(let [;; Serializtion format, must use same val for client + server:
      packer :edn ; Default packer, a good choice in most cases
      ;; (sente-transit/get-transit-packer) ; Needs Transit dep

      chsk-server
      (sente/make-channel-socket-server!
       (get-sch-adapter) {:packer packer :csrf-token-fn nil})

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

(defn wrap-formats [handler]
  (wrap-restful-format
   handler
   {:formats [:json-kw :transit-json :transit-msgpack]}))

(defroutes ring-routes
  (GET  "/"      ring-req (home-page ring-req))
  (GET  "/docs"  ring-req (fn [_]
                            (-> (response/ok (-> "docs/docs.md" io/resource slurp))
                                (response/header "Content-Type" "text/plain; charset=utf-8"))))
  (GET  "/chsk"  ring-req (ring-ajax-get-or-ws-handshake ring-req))
  (wrap-formats (POST "/json"  ring-req (fn [request] (response/ok
                                                       {:result (-> request :params)}))))
  (POST "/chsk"  ring-req (ring-ajax-post                ring-req))
  (POST "/login" ring-req (login-handler                 ring-req))
  (route/resources "/") ; Static files, notably public/main.js (our cljs target)
  (route/not-found "<h1>Page not found</h1>"))

;; (defroutes data


(def main-ring-handler
  "**NB**: Sente requires the Ring `wrap-params` + `wrap-keyword-params`
  middleware to work. These are included with
  `ring.middleware.defaults/wrap-defaults` - but you'll need to ensure
  that they're included yourself if you're not using `wrap-defaults`.
  You're also STRONGLY recommended to use `ring.middleware.anti-forgery`
  or something similar."
  (ring.middleware.defaults/wrap-defaults
   (-> ring-routes var wrap-reload wrap-webjars)
   (assoc-in ring.middleware.defaults/site-defaults [:security :anti-forgery] false)))
   ;; ring.middleware.defaults/site-defaults))

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
(defmethod -event-msg-handler :example/button1
  [{:as ev-msg :keys [?reply-fn]}]
  (println "button1 pushed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))

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
(defn start! [] (start-router!) (start-web-server! 3000) (start-example-broadcaster!))

(defn -main "For `lein run`, etc." [] (start!))

;; (println (str "nrepl" (env :nrepl-port)))

(comment
  (start!)
  (test-fast-server>user-pushes))
