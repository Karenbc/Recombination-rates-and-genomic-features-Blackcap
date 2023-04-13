BEGIN{
	OFS="\t"
}
{
	W=$2-$1
	r=$3
	cM=100*W*r
	cMb=r*1e+8
	Map+=100*W*r
	print $1,$2,cMb,Map,W
}

