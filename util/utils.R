performORA = function(set, all) {
	set = as.logical(set)
	all = as.logical(all)

	a = sum(set)
	b = sum(all) - a
	c = sum(!set)
	d = sum(!all) - c
	fisher.test(matrix(c(a, c, b, d),2,2), alternative='greater')$p.value
}
