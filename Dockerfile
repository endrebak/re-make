FROM openjdk:8-alpine

COPY target/uberjar/re-make.jar /re-make/app.jar

EXPOSE 3000

CMD ["java", "-jar", "/re-make/app.jar"]
