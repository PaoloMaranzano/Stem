# Archived: the Po Valley PM2.5 application in the paper

On 2026-09-27 every reference to the air quality application was removed from
the SC-STEM paper, because the application is going to change. The text as it
stood just before the removal is kept here in full:

| file | content |
|---|---|
| `main-before-removal.tex` | the manuscript, with Section 4 "Application: fine particulate matter over the Po Valley", the air quality motivation of the introduction, and the numbers of the application (d = 36, T = 1826) used in the tuning section |
| `supplement-before-removal.tex` | the supplement, with the section "Additional material for the application" and the numbers of the application in Section S2 |

Search the files for "Po Valley", "PM", "pollut", "episode" or "1{,}826" to
find the passages. The application script itself, `run-application.R`, is still
in `dev/replication` and runs on the `povalley` data shipped with the package.
