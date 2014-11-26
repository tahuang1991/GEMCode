

'''A function returning a string and using a local variable'''

def lastFirst(firstName, lastName):
        separator = ', '
	result = lastName + separator + firstName
        return result
str1 = lastFirst('A','B')
print "str1 ", str1, " type ",type(str1)
print "function return type ",type(lastFirst('A','B'))

print(lastFirst('Benjamin', 'Franklin'))
print(lastFirst('Andrew', 'Harrington'))

