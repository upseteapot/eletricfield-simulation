from matplotlib import pyplot as plt
import numpy as np


def sim_step(prev, cols, rows):
    new = np.zeros((rows, cols))
    for y in range(1, cols-2):
        for x in range(1, cols-2):
          new[y,x] += 1/4 * (prev[y-1,x] + prev[y,x-1] + prev[y+1,x] + prev[y,x+1])
    return new


if __name__ == "__main__":
    world_width  = 10.0 # in meters
    world_height = 10.0 # in meters
    world_res    = 0.05  # dx size
    
    cols = int(world_width/world_res)
    rows = int(world_height/world_res)
    
    world_cells  = np.zeros((rows, cols))
    
    for _ in range(100):
        for i in range(30,cols-30):
            world_cells[40,i] = 10
            world_cells[50,i] = -10
        world_cells = sim_step(world_cells, cols, rows)

    x_side = np.linspace(0, world_width, cols)
    y_side = np.linspace(0, world_height, rows)
    X, Y = np.meshgrid(x_side, y_side)
    
    plt.pcolormesh(X, Y, world_cells)
    plt.show()

