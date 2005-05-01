assign.groups <- function(object, groups)
{
  if(length(object$groups))
	  stop("tentative groups are already assigned.")

	n <- object$n
	
	n.groups <- length(groups)
	if(is.null(names(groups)))
	  group.names <- paste("Group", 1:n.groups)
	else
	  group.names <- names(groups)

	unassigned <- setdiff(1:n, unlist(groups))
	
	object$groups <- groups
	object$group.names <- group.names
	object$unassigned <- unassigned
	
	object
}


