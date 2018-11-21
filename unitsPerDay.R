# unitsPerDay.R
#
# calculate average number of ABC learning units to complete, per remaining
# day in the term.

nUnits <- 115     # Number of learning units on the map
nMilestones <- 3  # Number of "milestone" units - nothing to be done here

nDone <- 23        # Number of units I have completed by now (Nov11: 5-7 today!)

today <- Sys.Date()
lastDay <- as.Date("2018-12-06",'%Y-%m-%d')

nDays <- as.integer(lastDay - today)
avUnits <- (nUnits - nMilestones - nDone) / nDays

cat(sprintf("%1.1f units to complete per day in the remaining %d days.\n"
            , avUnits
            , nDays))

# [END]
