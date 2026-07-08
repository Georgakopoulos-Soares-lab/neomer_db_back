# Railway Deployment

This repository is ready for Railway deployment through the GitHub repository using the included Dockerfile.

## Required Railway variables

Set `NEOMERS_DUCK_DB_FILE` to the absolute path of the DuckDB database file inside the deployed container.

If you use a Railway volume, mount it at `/data` and place the database at:

```text
/data/staging.neomers.ddb
```

That path is already configured as the Docker default and can be overridden in Railway.

## Health check

Railway is configured to use:

```text
/healthcheck
```

The application must be able to open the DuckDB file during startup before the health check can pass.
