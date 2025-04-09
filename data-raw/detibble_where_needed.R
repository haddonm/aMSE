

library(tibble)

# Creating a tibble
my_tibble <- tibble(
  `year` = 1:5,
  `catch` = letters[1:5],
  `cpue` = c(TRUE, FALSE, TRUE, FALSE, TRUE)
)

# Printing the tibble
print(my_tibble)

# Accessing a column with a non-syntactic name
my_tibble$year
#Summary
#Tibbles are a powerful and user-friendly alternative to traditional 
#data frames, offering improved printing, subsetting, and error handling, 
#along with better integration with modern data manipulation packages like 
#dplyr. They make working with data in R more efficient and less error-prone.

