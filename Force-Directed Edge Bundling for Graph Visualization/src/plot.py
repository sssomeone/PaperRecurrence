#C:\\Users\\pengfei\\Desktop\\result.txt
import matplotlib.pyplot as plt

resultPath="C:\\Users\\pengfei\\Desktop\\result.txt"
with open(resultPath,"r") as file:
    while 1:
        line = file.readline()
        if not line:
            break
        xy=line[:-1].split(" ")
        plt.plot([float(xy[0]),float(xy[1])],[float(xy[2]),float(xy[3])])

plt.savefig('result.jpg')
plt.show()