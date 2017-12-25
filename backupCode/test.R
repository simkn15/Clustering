mylist<-list()
# enqueue
mylist <- c(mylist, list(1:5))
mylist <- c(mylist, list(5:10))
#dequeue
first <- mylist[[1]]
mylist <- mylist[-1]
