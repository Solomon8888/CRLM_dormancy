#!/usr/bin/env zsh

# zsh脚本共用终端运行函数
#
# 设计目标：
# - 构建、编译、清理这类固定流程只显示单行动态状态，不刷屏；
# - 真实命令输出写入temporary日志，失败时再打印，便于排查；
# - 成功时给出简洁的耗时统计，保持和R脚本的运行风格一致。

format_elapsed_seconds() {
  local total_seconds="${1:-0}"
  local minutes=$(( total_seconds / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02d:%02d" "$minutes" "$seconds"
}

run_with_spinner() {
  local label="$1"
  local log_file="$2"
  shift 2

  mkdir -p "$(dirname "$log_file")"
  : > "$log_file"

  local start_time
  start_time=$(date +%s)

  if [[ ! -t 1 ]]; then
    printf "… %s\n" "$label"
    "$@" >"$log_file" 2>&1
    local exit_status=$?
    local end_time elapsed
    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))

    if [[ "$exit_status" -eq 0 ]]; then
      printf "✓ %s | done in %s\n" "$label" "$(format_elapsed_seconds "$elapsed")"
      return 0
    fi

    printf "✗ %s | failed after %s\n" "$label" "$(format_elapsed_seconds "$elapsed")"
    printf "Log: %s\n" "$log_file"
    printf -- "---------------- command log ----------------\n"
    sed -n '1,220p' "$log_file"
    printf -- "---------------------------------------------\n"
    return "$exit_status"
  fi

  "$@" >"$log_file" 2>&1 &
  local command_pid=$!
  local frames=("⠋" "⠙" "⠹" "⠸" "⠼" "⠴" "⠦" "⠧" "⠇" "⠏")
  local frame_index=0
  local now elapsed end_time

  while kill -0 "$command_pid" 2>/dev/null; do
    now=$(date +%s)
    elapsed=$(( now - start_time ))
    printf "\r\033[2K%s %s | elapsed %s" \
      "${frames[$(( frame_index % ${#frames[@]} + 1 ))]}" \
      "$label" \
      "$(format_elapsed_seconds "$elapsed")"
    frame_index=$(( frame_index + 1 ))
    sleep 0.18
  done

  wait "$command_pid"
  local exit_status=$?
  end_time=$(date +%s)
  elapsed=$(( end_time - start_time ))

  if [[ "$exit_status" -eq 0 ]]; then
    printf "\r\033[2K✓ %s | done in %s\n" "$label" "$(format_elapsed_seconds "$elapsed")"
    return 0
  fi

  printf "\r\033[2K✗ %s | failed after %s\n" "$label" "$(format_elapsed_seconds "$elapsed")"
  printf "Log: %s\n" "$log_file"
  printf -- "---------------- command log ----------------\n"
  sed -n '1,220p' "$log_file"
  printf -- "---------------------------------------------\n"
  return "$exit_status"
}

print_step_header() {
  local title="$1"
  printf "\n%s\n" "$title"
  printf "%*s\n" "${#title}" "" | tr " " "-"
}
