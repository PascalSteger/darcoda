import hashlib
m = hashlib.md5()
#m.update("Nobody inspects")
#m.update(" the spammish repetition")
m.update(str(270*1e8+4))
print m.hexdigest()
