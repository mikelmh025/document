x = 10
def changeme():
    global x
    x = "works"
print x
changeme()
print x