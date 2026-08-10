#!/usr/bin/env bash
# Generic background job runner for GenoCore pipelines.
#
# Usage:
#   ./run_bg.sh start "your command here"     # Start a background job
#   ./run_bg.sh status                        # Check job status
#   ./run_bg.sh logs                          # Tail the log file
#   ./run_bg.sh stop                          # Stop the running job
#   ./run_bg.sh restart "your command here"   # Stop + restart with new command
#
# Examples:
#   ./run_bg.sh start "Rscript Genocore/run_genocore.R mydata.csv -cv 99 -o mydata"
#   ./run_bg.sh start "plink --vcf data.vcf.gz --make-bed --out data_qc"
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
LOG_DIR="$SCRIPT_DIR/logs"
PID_FILE="$SCRIPT_DIR/.daemon.pid"
LOG_NAME_FILE="$SCRIPT_DIR/.daemon.logname"

# ========== Helper functions ==========
get_log_file() {
  if [[ -f "$LOG_NAME_FILE" ]]; then
    cat "$LOG_NAME_FILE"
  else
    echo ""
  fi
}

is_running() {
  if [[ -f "$PID_FILE" ]]; then
    local pid=$(cat "$PID_FILE")
    if ps -p "$pid" >/dev/null 2>&1; then
      return 0
    else
      rm -f "$PID_FILE" "$LOG_NAME_FILE"
    fi
  fi
  return 1
}

# ========== Commands ==========
start() {
  if [[ $# -eq 0 ]]; then
    echo "Usage: $0 start \"command to run\""
    exit 1
  fi
  local cmd="$1"
  if is_running; then
    echo "Error: job already running (PID: $(cat "$PID_FILE"))"
    exit 0
  fi
  mkdir -p "$LOG_DIR"
  TEMP_LOG_FILE="${LOG_DIR}/run_$(date +%Y%m%d_%H%M%S).log"
  nohup bash -c "$cmd" >> "$TEMP_LOG_FILE" 2>&1 &
  local pid=$!
  echo "$pid" > "$PID_FILE"
  echo "$TEMP_LOG_FILE" > "$LOG_NAME_FILE"
  echo "Job started (PID: $pid)"
  echo "Log: $TEMP_LOG_FILE"
}

status() {
  local actual_log_file=$(get_log_file)
  if is_running; then
    echo "Running (PID: $(cat "$PID_FILE"))"
    echo "--- Last 20 log lines ---"
    if [[ -n "$actual_log_file" && -f "$actual_log_file" ]]; then
      tail -n 20 "$actual_log_file" || true
    else
      echo "Warning: log file not found"
    fi
  else
    echo "No job running"
    if [[ -n "$actual_log_file" && -f "$actual_log_file" ]]; then
      echo "--- Last run log tail ---"
      tail -n 10 "$actual_log_file" || true
    fi
  fi
}

stop() {
  if ! is_running; then
    echo "No job running"
    return 0
  fi
  local pid=$(cat "$PID_FILE")
  echo "Stopping job (PID: $pid)..."
  kill -TERM "$pid" 2>/dev/null || true
  local count=0
  while ps -p "$pid" >/dev/null 2>&1 && [[ $count -lt 30 ]]; do
    sleep 1
    count=$((count + 1))
  done
  if ps -p "$pid" >/dev/null 2>&1; then
    echo "Force killing (PID: $pid)"
    kill -KILL "$pid" 2>/dev/null || true
  fi
  rm -f "$PID_FILE" "$LOG_NAME_FILE"
  echo "Job stopped"
}

restart() {
  if [[ $# -eq 0 ]]; then
    echo "Usage: $0 restart \"command to run\""
    exit 1
  fi
  local cmd="$1"
  echo "Restarting..."
  stop || true
  sleep 2
  start "$cmd"
}

logs() {
  local actual_log_file=$(get_log_file)
  mkdir -p "$LOG_DIR"
  if [[ -n "$actual_log_file" && -f "$actual_log_file" ]]; then
    echo "Tailing log (Ctrl+C to exit):"
    tail -n 200 -f "$actual_log_file"
  else
    echo "Warning: no log file found"
  fi
}

# ========== Main ==========
cmd="${1:-start}"
case "$cmd" in
  start)
    shift
    start "$@"
    ;;
  status)
    status
    ;;
  stop)
    stop
    ;;
  restart)
    shift
    restart "$@"
    ;;
  logs)
    logs
    ;;
  *)
    echo "Usage: $0 {start|status|stop|restart|logs} [command]"
    exit 1
    ;;
esac
