student <- setRefClass("student", fields = list(name = "character", age = "numeric", GPA = "numeric"))
# ^generator func       ^className  ^classFields
# now student() is our generator function which can be used to create new objects
s <- student(name = "John", age = 21, GPA = 3.5)
#Only a single copy exist and all variables reference to the same copy. Hence the name, reference.
b <- s
b$name <- "simon"
s
b
# create reference object a and assign a’s copy to b
a <- student(name = "John", age = 21, GPA = 3.5)
b <- a$copy()
b$name <- "Paul"
a
b

you <- setRefClass("test", fields = list(data = "list"))
xx <- you(data = list(1,2,3))
xx$getClass()
