##############################################
# Sliding-window π and Tajima's D (SR vs ST)
##############################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Paths to input files (thetasWindow *.pestPG). Gzipped OK.
file_sr <- "filename.thetasWindow.gz.pestPG"
file_st <- "filename.thetasWindow.gz.pestPG"

# Windowing parameters
WINDOW_BP <- 500000   # 500 kb window
STEP_BP   <-  50000   # 50 kb step

# Chromosome order & labels (facet order + human names)
chrom_order <- c("NC_046683.1", "NC_046679.1", "NC_046680.1", "NC_046681.1")
chrom_labels <- c(
  "NC_046683.1" = "X chromosome",
  "NC_046679.1" = "2nd chromosome",
  "NC_046680.1" = "3rd chromosome",
  "NC_046681.1" = "4th chromosome"
)

# Mark SR-chromosome inversion breakpoints
zoom_chr   <- "NC_046683.1"
landmarks  <- c(44280095, 46899044, 47586979, 56998335, 64000000, 67500000)

# Read + standardize a single dataset to a minimal, consistent schema.
# Expects columns: Chr, WinCenter, tP, Tajima
read_thetas <- function(path, sample_label) {
  if (!file.exists(path)) {
    stop(sprintf("Input not found: %s", path))
  }
  # read.table handles gz; keep character then coerce explicitly
  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   colClasses = "character")
  
  # Coerce and select only what we truly need downstream
  df %>%
    transmute(
      sample     = sample_label,
      chromosome = Chr,
      coord      = suppressWarnings(as.numeric(WinCenter)),
      pi         = suppressWarnings(as.numeric(tP)) / 50000,
      tajimaD    = suppressWarnings(as.numeric(Tajima))
    )
}

# Compute sliding windows per (sample, chromosome)
sliding_window_summary <- function(df, window_bp, step_bp) {
  # Require numeric coords
  df <- df %>% filter(!is.na(coord))
  
  if (nrow(df) == 0) return(tibble())
  
  chr_min <- min(df$coord, na.rm = TRUE)
  chr_max <- max(df$coord, na.rm = TRUE)
  
  # Last start not past (max - window); if window > span, keep one window at chr_min
  seq_end <- max(chr_max - window_bp, chr_min)
  starts  <- seq(chr_min, seq_end, by = step_bp)
  
  windows <- lapply(starts, function(start) {
    end   <- start + window_bp
    slice <- df %>% filter(coord >= start, coord < end)
    if (nrow(slice) == 0) return(NULL)
    
    tibble(
      start        = start,
      end          = end,
      mid          = start + window_bp / 2,
      pi_mean      = mean(slice$pi, na.rm = TRUE),
      tajimaD_mean = mean(slice$tajimaD, na.rm = TRUE),
      chromosome   = dplyr::first(slice$chromosome),
      sample       = dplyr::first(slice$sample)
    )
  })
  
  bind_rows(windows)
}

base_theme <- theme_minimal(base_size = 14) +
  theme(
    strip.text    = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 12),
    axis.title    = element_text(size = 14),
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(0.35, "lines")
  )

# ---- 3) Load & combine data ----
data_sr <- read_thetas(file_sr, "SR")
data_st <- read_thetas(file_st, "ST")

data_all <- bind_rows(data_sr, data_st) %>%
  mutate(chromosome = factor(chromosome, levels = chrom_order))

# ---- 4) Sliding windows by (sample, chromosome) ----
windowed <- data_all %>%
  group_by(sample, chromosome, .drop = FALSE) %>%
  group_split() %>%
  lapply(sliding_window_summary, window_bp = WINDOW_BP, step_bp = STEP_BP) %>%
  bind_rows() %>%
  mutate(chromosome = factor(chromosome, levels = chrom_order))

# ---- 5) π plot (sliding window) ----
plot_pi <- ggplot(windowed, aes(x = mid, y = pi_mean, color = sample)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(SR = "red", ST = "blue")) +
  facet_grid(. ~ chromosome, scales = "free_x", space = "free_x",
             labeller = labeller(chromosome = chrom_labels)) +
  labs(x = NULL, y = expression(pi)) +
  coord_cartesian(ylim = c(0, 0.025)) +
  base_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- 6) Tajima's D plot (sliding window) ----
plot_td <- ggplot(windowed, aes(x = mid, y = tajimaD_mean, color = sample)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(SR = "red", ST = "blue")) +
  facet_grid(. ~ chromosome, scales = "free_x", space = "free_x",
             labeller = labeller(chromosome = chrom_labels)) +
  labs(x = NULL, y = "Tajima's D") +
  coord_cartesian(ylim = c(-2, 0.1)) +
  base_theme +
  theme(
    strip.text   = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- 7) SR/ST π ratio plot ----
pi_ratio <- windowed %>%
  group_by(chromosome, mid, sample) %>%
  summarise(pi_mean = mean(pi_mean, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sample, values_from = pi_mean) %>%
  mutate(
    pi_ratio  = ifelse(!is.na(ST) & ST > 0, SR / ST, NA_real_),
    chromosome = factor(chromosome, levels = chrom_order)
  )

plot_ratio <- ggplot(pi_ratio, aes(x = mid, y = pi_ratio)) +
  geom_line(linewidth = 0.9, color = "black") +
  facet_grid(. ~ chromosome, scales = "free_x", space = "free_x",
             labeller = labeller(chromosome = chrom_labels)) +
  labs(x = "Coordinate", y = expression(SR/ST~pi)) +
  base_theme +
  theme(strip.text = element_blank())

# ---- 8) Combine panels (stacked) ----
figure_main <- plot_pi / plot_td / plot_ratio + plot_layout(heights = c(1, 1, 1))

# Print to device
print(figure_main)

# ---- 9) Zoomed view on X chromosome----
win_zoom   <- filter(windowed, chromosome == zoom_chr)
ratio_zoom <- filter(pi_ratio, chromosome == zoom_chr)

add_landmarks <- function(p) {
  if (length(landmarks) == 0) return(p)
  p + geom_vline(xintercept = landmarks, linewidth = 0.5)
}

plot_pi_zoom <- ggplot(win_zoom, aes(x = mid, y = pi_mean, color = sample)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(SR = "red", ST = "blue")) +
  labs(x = NULL, y = expression(pi)) +
  coord_cartesian(ylim = c(0, 0.025)) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) |>
  add_landmarks()

plot_td_zoom <- ggplot(win_zoom, aes(x = mid, y = tajimaD_mean, color = sample)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(SR = "red", ST = "blue")) +
  labs(x = NULL, y = "Tajima's D") +
  coord_cartesian(ylim = c(-2, 0.25)) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) |>
  add_landmarks()

plot_ratio_zoom <- ggplot(ratio_zoom, aes(x = mid, y = pi_ratio)) +
  geom_line(linewidth = 0.9, color = "black") +
  labs(x = "Coordinate", y = expression(SR/ST~pi)) +
  base_theme |>
  add_landmarks()

figure_zoom <- plot_pi_zoom / plot_td_zoom / plot_ratio_zoom +
  plot_layout(heights = c(1, 1, 1))

print(figure_zoom)
