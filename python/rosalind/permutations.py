#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Input: number n<=7
Output:
	permutation counts
	perm1
	perm2
	...

2017.6.21 by xnm
'''
f = open("rosalind_perm.txt","r")
n = int(f.readline())
f.close()

all_perm = []
num = range(1,n+1)
def permute(numList, low=0):
	if low + 1 >= len(numList):
		all_perm.append(tuple(numList))
	else:
		permute(numList, low + 1)
		for i in range(low + 1, len(numList)):
			numList[low],numList[i] = numList[i], numList[low]
			permute(numList, low+1)
		numList[low],numList[i] = numList[i], numList[low]

permute(num)


ct = len(all_perm)
print ct
for i in all_perm:
	for j in i:
		print j,
	print "hello"

