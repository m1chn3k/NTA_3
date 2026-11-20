FROM gcc:14
WORKDIR /app
COPY NTA_lab3/*.cpp NTA_lab3/*.h ./
RUN g++ -std=c++17 -O2 -pthread *.cpp -o lab3
ENTRYPOINT ["./lab3"]
