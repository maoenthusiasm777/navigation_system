import matplotlib.pyplot as plt

with open("./zuhe_pos_vel_yaw.txt") as f:
    datas = []
    for line in f.readlines():
        data = []
        for pos_vel in line.split():
            data.append(float(pos_vel))
        datas.append(data)
plt.plot([pos_vel[1] for pos_vel in datas],[pos_vel[0] for pos_vel in datas])
plt.show()