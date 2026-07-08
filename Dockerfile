FROM golang:1.23-bookworm AS builder

WORKDIR /app

COPY go.mod go.sum ./
RUN go mod download

COPY server.go ./
RUN CGO_ENABLED=1 GOOS=linux GOARCH=amd64 go build -trimpath -ldflags="-s -w" -o /out/neomer_server ./server.go

FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates libgomp1 libstdc++6 socat \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /out/neomer_server ./neomer_server

ENV GIN_MODE=release \
    NEOMERS_DUCK_DB_FILE=/data/staging.neomers.ddb

EXPOSE 8080

CMD ["sh", "-c", "set -eu; export PORT=\"${PORT:-8080}\"; if [ \"$PORT\" = \"8080\" ]; then exec ./neomer_server; else socat TCP-LISTEN:${PORT},fork,reuseaddr TCP:127.0.0.1:8080 & exec ./neomer_server; fi"]
