//! Centralized styling for CLI output (palette adapted from escapepod-rs).
//!
//! Every helper returns an already-rendered `String` — ANSI when color is
//! enabled, plain text otherwise. Color is suppressed automatically when:
//!
//! - `NO_COLOR` is set (<https://no-color.org/>)
//! - `CLICOLOR=0` is set
//! - stderr is not a TTY (e.g. piped or redirected)
//!
//! `CLICOLOR_FORCE=1` overrides the TTY check and always emits ANSI.
//!
//! Status output goes to stderr, so the gate keys off *stderr*'s terminal
//! status; machine-readable payloads (TSV/JSON/Stockholm) are written
//! unstyled to stdout regardless.

// Shared theme: a few helpers are kept for a consistent vocabulary even when not
// every command uses each one yet.
#![allow(dead_code)]

use owo_colors::OwoColorize;
use std::fmt::Display;
use std::io::IsTerminal;
use std::sync::OnceLock;

static USE_COLOR: OnceLock<bool> = OnceLock::new();

/// Whether ANSI styling should be emitted (cached on first call).
pub fn use_color() -> bool {
    *USE_COLOR.get_or_init(|| {
        // Force-on wins over everything.
        if matches!(
            std::env::var("CLICOLOR_FORCE").as_deref(),
            Ok(v) if v != "0" && !v.is_empty()
        ) {
            return true;
        }
        // NO_COLOR: any value, including empty, disables color.
        if std::env::var_os("NO_COLOR").is_some() {
            return false;
        }
        // CLICOLOR=0 disables.
        if matches!(std::env::var("CLICOLOR").as_deref(), Ok("0")) {
            return false;
        }
        // Status prints go to stderr; gate on stderr's terminal status.
        std::io::stderr().is_terminal()
    })
}

/// Action verbs like "Scanning", "Searching", "Aligning", "Found".
pub fn action<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.green().bold())
    } else {
        s.to_string()
    }
}

/// File paths.
pub fn path<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.cyan())
    } else {
        s.to_string()
    }
}

/// Important counts/numbers (hits found, models scanned, …).
pub fn count<T: Display>(n: T) -> String {
    if use_color() {
        format!("{}", n.green())
    } else {
        n.to_string()
    }
}

/// Percentages.
pub fn percentage<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.cyan())
    } else {
        s.to_string()
    }
}

/// Inline labels like "Output:", "Threshold:".
pub fn label<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.bold())
    } else {
        s.to_string()
    }
}

/// Section headers ("Modification database").
pub fn header<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.bold().cyan())
    } else {
        s.to_string()
    }
}

/// Key names in key/value pairs.
pub fn key<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.blue())
    } else {
        s.to_string()
    }
}

/// Values in key/value pairs.
pub fn value<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.cyan())
    } else {
        s.to_string()
    }
}

/// Dimmed secondary text / units.
pub fn dim<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.dimmed())
    } else {
        s.to_string()
    }
}

/// Warning label.
pub fn warning_label<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.yellow().bold())
    } else {
        s.to_string()
    }
}

/// Warning messages/values.
pub fn warning<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.yellow())
    } else {
        s.to_string()
    }
}

/// Error messages/values.
pub fn error<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.red())
    } else {
        s.to_string()
    }
}

/// Note / "Skipped" prefixes.
pub fn note_label<T: Display>(s: T) -> String {
    if use_color() {
        format!("{}", s.yellow())
    } else {
        s.to_string()
    }
}
