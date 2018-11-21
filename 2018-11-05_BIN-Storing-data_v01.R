# ==   1.3  Task: submit for credit (part 1/2)  ================================


#    Write and submit code that adds another philosopher to the datamodel:
#       Immanuel Kant, (1724 - 1804), Enlightenment Philosophy.
#       Works: Critique of Pure Reason (1781), Critique of Judgement (1790)

(newPhil <- data.frame(
  id = autoincrement(philDB$person),
  name = "Immanuel Kant",
  born = "1724",
  died = "1804",
  school = "Enlightenment Philosophy",
  stringsAsFactors = FALSE))

philDB$person <- rbind(philDB$person, newPhil)

(newBooks <- data.frame(
  id = c(autoincrement(philDB$books),autoincrement(philDB$books)+1),
  title = c("Critique of Pure Reason", "Critique of Judgement"),
  published = c("1781", "1790"),
  stringsAsFactors = FALSE))

philDB$books <- rbind(philDB$books, newBooks)

(newWorks <- data.frame(
  id = c(autoincrement(philDB$works), autoincrement(philDB$works)+1),
  personID = c(grep("Immanuel Kant", philDB$person$name),
               grep("Immanuel Kant", philDB$person$name)),
  bookID = c(grep("Critique of Pure Reason", philDB$books$title),
             grep("Critique of Judgement", philDB$books$title)),
  stringsAsFactors = FALSE))

philDB$works <- rbind(philDB$works, newWorks)

#    Write and submit code that lists the books in alphabetical order,
#    followed by the author and the year of publishing. Format your output like:
#    "Analects" - Kongzi (220 BCE)
#    Show the result.

alphabeticBookID <- order(philDB$books$title)
for (ID in alphabeticBookID) {
  workID <- which(philDB$works$bookID == ID)
  philID <- philDB$works$personID[workID]

  cat(sprintf("\"%s\" - ", philDB$books$title[ID]))
  cat(sprintf("%s (%s) \n",
              philDB$person$name[philID],
              philDB$books$published[ID]))
}
