import pickle

# need to deal with Python 2 objects explicitly in Py3.
def pickleRead(handle):
    read= pickle.load(handle, encoding='latin-1')

    # need to covert the strings to be in the right format.
    for i in range(len(read['chipNames'])):
        temp= read['chipNames'][i]
        temp2= []
        for index, value in enumerate(temp):
            temp2.append(value.decode("utf-8"))

        #print(temp2)
        read['chipNames'][i]= temp2
    return read
