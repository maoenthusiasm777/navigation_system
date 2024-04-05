import matplotlib.pyplot as plt

with open("./zuhe_pos_vel_yaw1.txt") as f:
    datas = []
    for line in f.readlines():
        data = []
        for pos_vel in line.split():
            data.append(float(pos_vel))
        datas.append(data)

with open("./gps_data.txt") as f:
    gnssdatas = []
    for line in f.readlines():
        gnssdata = []
        for pos_vel in line.split():
            gnssdata.append(float(pos_vel))
        gnssdatas.append(gnssdata)

plt.plot([pos_vel[9] for pos_vel in gnssdatas],[pos_vel[7] for pos_vel in gnssdatas])
plt.plot([pos_vel[1] for pos_vel in datas],[pos_vel[0] for pos_vel in datas])
plt.show()