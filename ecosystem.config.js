module.exports = {
  apps: [
    {
      name: "sybr-api",
      script: "sybr_api.py",
      interpreter: "/home/pranjal.p/.conda/envs/sybr/bin/python",
      args: "--port 8001",
      cwd: "/home/pranjal.p/sybr/Sybr/",
      instances: 1,
      autorestart: true,
      watch: false,
      max_restarts: 5,
      restart_delay: 3000,
      env: {
        SYBR_PIPELINE_DIR: "/home/pranjal.p/sybr/Sybr",
        SYBR_JOBS_DIR: "/home/pranjal.p/sybr/Sybr/jobs",
        SYBR_DB_PATH: "/home/pranjal.p/sybr/Sybr/api/sybr_api.db",
        SYBR_API_HOST: "0.0.0.0",
        SYBR_API_PORT: "8001",
        SYBR_MAX_JOBS: "2",
      },
      log_date_format: "YYYY-MM-DD HH:mm:ss",
      out_file: "/home/pranjal.p/.pm2/logs/sybr-api-out.log",
      error_file: "/home/pranjal.p/.pm2/logs/sybr-api-error.log",
      merge_logs: false,
    },
  ],
};
