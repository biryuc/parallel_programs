
import os
import numpy as np
import matplotlib.pyplot as plt
import csv

def test_calc_pi_iter(N):

    iter_arr = []
    tmp = 10
    for i in range(N):
        tmp = tmp * 10
        iter_arr.append(tmp)
    for i in range(N):
        str_ = str(float(iter_arr[i])) + " " + str(int(8))
        cmd = "C:\\Users\\Никита\\source\\repos\\calculate_pi\\x64\\Release\\calculate_pi.exe" + " " + str_ + " " + "1>>test_calc_pi.csv"
        os.system(cmd)
        print(i)

def test_calc_pi_thread():
    n = 0
    for i in range(8):
        n = n+1
        str_ = str(float(100000000)) + " " + str(int(n))
        cmd = "C:\\Users\\Никита\\source\\repos\\calculate_pi\\x64\\Release\\calculate_pi.exe" + " " + str_ + " " + "1>>test_calc_pi_thread.csv"
        os.system(cmd)
        print(n)

if __name__ == '__main__':
    test = 0
    N = 5

    test_calc_pi_iter(N)
    test_calc_pi_thread()


    filenames   = ["test_calc_pi.csv","test_calc_pi_thread.csv"]
    iter = []
    time = []
    time_omp = []
    res = []
    num = []

    iter2 = []
    time2 = []
    time_omp2 = []
    res2 = []
    num2 = []

    boost = []
    effect = []
    counter = 0

    with open(filenames[0], 'r',newline="\n") as csvfile:
        lines = csv.reader(csvfile, delimiter='\t')
        for row in lines:
            iter.append(int(row[0]))
            time.append(float(row[1]))
            time_omp.append((float(row[2])))
            res.append(int(row[3]))
            num.append((int(row[4])))

    csvfile.close()
    f = open(filenames[0], 'w+')
    f.seek(0)
    f.close()

    f = open("res.csv",'w')
    with open(filenames[1], 'r',newline="\n") as csvfile:
        lines = csv.reader(csvfile, delimiter='\t')
        for row in lines:
            iter2.append(int(row[0]))
            f.write(str(row[0]))
            f.write('\t')
            time2.append(float(row[1]))
            f.write(str(row[1]))
            f.write('\t')
            time_omp2.append((float(row[2])))
            f.write(str(row[2]))
            f.write('\t')
            res2.append(int(row[3]))
            f.write(str(row[3]))
            f.write('\t')
            num2.append((int(row[4])))
            f.write(str(row[4]))
            f.write('\n')
    csvfile.close()
    f.close()
    f = open(filenames[1], 'w+')
    f.seek(0)
    f.close()

    average = sum(time2) / len(time2)


    for i in range(len(time2)):
        boost.append(average/time_omp2[i])

    for i in range(len(time2)):
        effect.append(boost[i]*100/(i+1))

    plt.figure(figsize=(7, 7))
    plt.subplot(2, 2, 1)
    plt.plot(num2, time_omp2, color='b', label="time omp")
    plt.xlabel('treads')
    plt.ylabel('time')
    plt.title('Work time', fontsize=10)
    plt.subplot(2, 2, 2)
    plt.plot(iter, time, color='g', label="time")
    plt.plot(iter, time_omp, color='r',label="time omp")
    plt.xlabel('iter')
    plt.ylabel('time')
    plt.title('Work time', fontsize=10)
    plt.subplot(2, 2, 3)
    plt.plot(num2, boost, color='g', label="boost")
    plt.xlabel('treads')
    plt.ylabel('time')
    plt.title('Ускорение', fontsize=10)
    plt.subplot(2, 2, 4)
    plt.plot(num2, effect, color='r', label="effect")
    plt.xlabel('treads')
    plt.ylabel('effect')
    plt.title('Эффективность', fontsize=10)
    plt.grid()
    plt.legend()
    plt.show()





