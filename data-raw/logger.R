# ==============================================================================
# data-raw/logger.R
#
# Lightweight tee-logger: mirrors cat() / message() / warning() output to
# both the console and a timestamped log file simultaneously.
#
# USAGE:
#   source("data-raw/logger.R")
#
#   log <- start_log("data-raw/logs")   # opens log file, returns path
#   # ... run your analysis ...
#   stop_log()                          # restores normal output
#
# What gets captured:
#   - cat() / print() output        via sink(split = TRUE, type = "output")
#   - message() output              via a globalCallingHandler
#   - warning() output              via a globalCallingHandler
#
# Note on R's message sink limitation:
#   R does not allow split = TRUE on the message connection
#   (it raises "cannot split the message connection"). Therefore messages and
#   warnings are captured by installing global calling handlers that write a
#   copy to the log file while letting the condition propagate normally so it
#   still appears on the console. stdout is handled with the normal split sink.
#
# If stop_log() is not called (e.g. the script errors out), call stop_log()
# manually or restart R to clear the sink.
# ==============================================================================


# ==============================================================================
# Internal state
# ==============================================================================
.log_env <- new.env(parent = emptyenv())
.log_env$active       <- FALSE
.log_env$log_path     <- NULL
.log_env$con          <- NULL
.log_env$old_handlers <- NULL   # saved globalCallingHandlers state


# ==============================================================================
# .write_to_log() — internal helper
# ==============================================================================

#' Write text to the active log file connection (silent on failure)
#' @noRd
.write_to_log <- function(text) {
  if (!.log_env$active || is.null(.log_env$con)) return(invisible(NULL))
  tryCatch(
    cat(text, file = .log_env$con, append = TRUE),
    error = function(e) NULL
  )
}


# ==============================================================================
# start_log()
# ==============================================================================

#' Open a tee-log file and start capturing console output
#'
#' Creates \code{log_dir} if it does not exist, opens a timestamped log file,
#' installs a split stdout sink, and registers global calling handlers for
#' messages and warnings so that every line of output is mirrored to the log.
#'
#' @param log_dir       Character. Directory for the log file.
#'   Default \code{"data-raw/logs"}.
#' @param timestamp_fmt Character. \code{format()} pattern for the filename.
#'   Default \code{"\%Y\%m\%d_\%H\%M\%S"}.
#'
#' @return The full path to the opened log file, invisibly.
start_log <- function(log_dir       = "data-raw/logs",
                      timestamp_fmt = "%Y%m%d_%H%M%S") {

  if (.log_env$active) {
    warning("A log is already active (", .log_env$log_path, "). ",
            "Call stop_log() before starting a new one.", call. = FALSE)
    return(invisible(.log_env$log_path))
  }

  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

  ts       <- format(Sys.time(), timestamp_fmt)
  log_path <- file.path(log_dir, paste0("update_network_", ts, ".log"))
  con      <- file(log_path, open = "wt")

  # ── stdout: split = TRUE mirrors cat()/print() to both console and file ────
  sink(con, append = FALSE, split = TRUE, type = "output")

  # ── messages & warnings: global calling handlers write a copy to the file ──
  # R does not support split = TRUE on the message connection, so we use
  # globalCallingHandlers (available since R 4.0) instead.
  # Each handler writes to the file then returns NULL, which lets the condition
  # propagate normally so it still appears on the console.
  .log_env$old_handlers <- globalCallingHandlers(

    message = function(m) {
      .write_to_log(conditionMessage(m))
      NULL   # let the default handler print to console
    },

    warning = function(w) {
      .write_to_log(paste0("Warning: ", conditionMessage(w), "\n"))
      NULL   # let the default handler print to console
    }
  )

  .log_env$active   <- TRUE
  .log_env$log_path <- log_path
  .log_env$con      <- con

  cat(sprintf("[log] Started : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("[log] File    : %s\n\n", log_path))

  invisible(log_path)
}


# ==============================================================================
# stop_log()
# ==============================================================================

#' Close the active log file and restore normal output
#'
#' Removes the split stdout sink and the global message/warning calling
#' handlers installed by \code{start_log()}. Calling this when no log is
#' active is a silent no-op.
#'
#' @return \code{NULL}, invisibly.
stop_log <- function() {

  if (!.log_env$active) return(invisible(NULL))

  cat(sprintf("\n[log] Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("[log] Saved to: %s\n",   .log_env$log_path))

  # Restore the previous global calling handlers
  tryCatch(
    globalCallingHandlers(.log_env$old_handlers),
    error = function(e) NULL
  )

  # Remove the stdout split sink
  tryCatch(sink(type = "output"), error = function(e) NULL)

  # Close the file connection
  tryCatch(close(.log_env$con), error = function(e) NULL)

  .log_env$active       <- FALSE
  .log_env$log_path     <- NULL
  .log_env$con          <- NULL
  .log_env$old_handlers <- NULL

  invisible(NULL)
}


# ==============================================================================
# log_path()
# ==============================================================================

#' Return the path of the currently active log file, or NULL
log_path <- function() .log_env$log_path
