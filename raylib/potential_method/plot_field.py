from matplotlib import pyplot as plt

X = list()
Y = list()

with open("field", "r") as file:
    for line in file:
        x_i, y_i = line.split(" ")
        X.append(float(x_i))
        Y.append(float(y_i))

plt.plot(X, Y)
plt.title("Electric field");
plt.show()

X = list()
Y = list()

with open("potential", "r") as file:
    for line in file:
        x_i, y_i = line.split(" ")
        X.append(float(x_i))
        Y.append(float(y_i))

plt.plot(X, Y)
plt.title("Electric potential");
plt.show()

