from sys import stderr
from numpy import *
import gurobipy as gp
from gurobipy import GRB
from gurobipy import read

examIDs = set({})
studIDs = set({})

# reading the slo instance
with open('C:/Users/asus/Desktop/ETP_problem/DODMproject2023/instances/test.slo') as file_slo:
    num_timeslots = tmax = int(file_slo.readline())

# reading the stu instance
with open('C:/Users/asus/Desktop/ETP_problem/DODMproject2023/instances/test.stu') as file_stu:
    stu_lines = file_stu.readlines()
    # print(f"{stu_lines=}", file=stderr)
    studs_enrolled_to_exam = {}
    for line in stu_lines:
        ID_stud, ID_exam = line.strip().split()
        ID_exam = int(ID_exam)
        examIDs.add(ID_exam)
        studIDs.add(ID_stud)
        if ID_exam in studs_enrolled_to_exam:
            if ID_stud in studs_enrolled_to_exam[ID_exam]:
                print(f"WARNING: double registration of student {ID_stud} to exam {ID_exam}!", file=stderr)
            else:
                studs_enrolled_to_exam[ID_exam].add(ID_stud)
        else:
            studs_enrolled_to_exam[ID_exam] = {ID_stud}
print(f"{examIDs=}", file=stderr)
print(f"{studIDs=}", file=stderr)
print(f"{studs_enrolled_to_exam=}", file=stderr)
num_exams = len(studs_enrolled_to_exam)

# reading the exm instance
with open('C:/Users/asus/Desktop/ETP_problem/DODMproject2023/instances/test.exm') as file_exm:
    exm_lines = file_exm.readlines()
    print(f"{exm_lines[:-1]=}", file=stderr)
    num_studs_enrolled_to_exam = {}
    for line in exm_lines[:-1]:
        ID_exam, num_enrolls = list(map(int, line.strip().split()))
        if ID_exam in num_studs_enrolled_to_exam:
            print(f"WARNING: double declaration of the number of enrollments to exam {ID_exam}!", file=stderr)
        num_studs_enrolled_to_exam[ID_exam] = num_enrolls
    print(f"{num_studs_enrolled_to_exam=}", file=stderr)

def check_consistency():
    for el in examIDs:
        if len(studs_enrolled_to_exam[el]) != num_studs_enrolled_to_exam[el]:
            print(
                f"Error: inconsistency detected! The number of students enrolled to exam {el} is {num_studs_enrolled_to_exam[el]}. However, the students enrolled according to the file are actually these:\n{studs_enrolled_to_exam[el]}")
            exit(1)

#check_consistency()

# creating parameters
# p = exponential contribute to the penalty
P = [1, 2, 3, 4, 5]
p = {}
for d in P:
    p[d] = 2 ** (5 - d)
print(p)

# E = examIDs
E = [x for x in range(1, 1 + num_exams)]
S = studIDs
T = [x for x in range(1, tmax + 1)]
print(E)
print(S)
print(T)
nexams = len(E)
nstudents = len(S)

# creating n[i, j]
n = {}
for i in examIDs:
    for j in examIDs:
        nstudconf = 0
        for stud in studIDs:
            if stud in studs_enrolled_to_exam[i] and stud in studs_enrolled_to_exam[j] and i != j:
                nstudconf += 1
        n[i, j] = nstudconf
# print(n)

# Setting up the environment
env = gp.Env(
    "C:/Users/asus/Desktop/ETP_problem/log/ETP_Problem_test.log")  # log file is what is printed in the console
env.start()

# Environment parameters
env.setParam("Threads", 1)  # Controls the number of threads to apply to parallel algorithms
env.setParam("Presolve", 1)  # Controls the presolve level. conservative (1).
env.setParam("MIPGap", 1e-4)
env.setParam('Method', 0)  # Algorithm used to solve the initial root relaxation of a MIP model. 0=primal simplex.
env.setParam("TimeLimit", 600)  # 10 minutes time limit
env.setParam("PreSparsify", 1) # to reduce the memory used
# env.setParam("MultiObjPre", 0)
# env.setParam("OutputFlag", 0)

# Initialization of a new model in the defined environment
model = gp.Model("ETP_Problem", env=env)

# VARIABLES
x = model.addVars(E, T, vtype=GRB.BINARY, name="x")  # equals 1 if exam i is schedule at time t, 0 otherwise
y = model.addVars(E, E, name="y")  # the penalty with exam i and j

# CONSTRAINTS

# every exam is assigned precisely once
model.addConstrs(gp.quicksum(x[i, t] for t in T) == 1 for i in E)
model.update()

# there can't be two conflicting exams in the same timeslot
model.addConstrs(x[i, t] + x[j, t] <= 1 for i in E for j in E for t in T if i < j and n[i, j] > 0)
model.update()

# constrain for defining the auxiliary variable
for d in P:
    for i in E:
        for j in E:
            for t in T:
                if i < j:
                    if t + d <= tmax and n[i, j] > 0:
                        model.addConstr(y[i, j] >= p[d] * x[i, t] + p[d] * x[j, t + d] - p[d])
                        model.addConstr(y[i, j] >= p[d] * x[i, t + d] + p[d] * x[j, t] - p[d])

# OBJECTIVE FUNCTION
#minimize the total penalty
model.setObjective(gp.quicksum(y[i, j] for i in E for j in E), GRB.MINIMIZE)

# Optimization
model.optimize()

print("\nObjective function value = {}".format(round(model.objVal, 5)))

# Reporting results
# print(model.Status)
# if model.Status == GRB.OPTIMAL:
#    print("\nOptimal objective function value = {}".format(round(model.objVal, 5)))
# print("{} = {}".format(x.VarName, x.getAttr("X")))
# print("{} = {}".format(y.VarName, y.getAttr("X")))
# else:
#    if model.Status == GRB.INFEASIBLE:
#        print("Infeasible")
#        model.computeIIS()
#        exit(0)
#    else:
#        print("Unbounded")

timeslot_of_exam = [None] * (num_exams + 1)
for i in E:
    for t in T:
        if bool(x[i, t].X):
            timeslot_of_exam[i] = t
            continue
print(f"timeslot_of_exam={timeslot_of_exam[1:]}", file=stderr)

exams_at_timeslot = [[] for _ in range(num_timeslots + 1)]
for ID_exam in range(1, 1 + num_exams):
    exams_at_timeslot[timeslot_of_exam[ID_exam]].append(ID_exam)
print(f"exams_at_timeslot={exams_at_timeslot[1:]}", file=stderr)


def check_feasible():
    for t in range(1, 1 + num_timeslots):
        activity = {s: None for s in studIDs}
        for e in exams_at_timeslot[t]:
            for s in studs_enrolled_to_exam[e]:
                if activity[s] is None:
                    activity[s] = e
                else:
                    print(
                        f"Ahi: the solution proposed is not feasible since at timeslot {t} the student {s} should attend both the exam {activity[s]} and the exam {e}!",
                        file=stderr)
                    return False
    return True


print(f"{check_feasible()=}")


def penality():
    p = 0.0
    students_active_in_timewind = [set({}), set({}), set({}), set({}), set({}), set({})]
    for t in range(1, 1 + num_timeslots):
        students_active_in_timewind[t % 6] = set({})
        activity = {s: None for s in studIDs}
        for e in exams_at_timeslot[t]:
            for s in studs_enrolled_to_exam[e]:
                if activity[s] is None:
                    activity[s] = e
                    students_active_in_timewind[t % 6].add(s)
                else:
                    print(
                        f"Ahi: the solution proposed is not feasible since at timeslot {t} the student {s} should attend both the exam {activity[s]} and the exam {e}!",
                        file=stderr)
                    exit(1)
            for time_shift in range(1, 6):
                p += 2 ** (5 - time_shift) * len(
                    studs_enrolled_to_exam[e].intersection(students_active_in_timewind[(t - time_shift) % 6])) / len(
                    studIDs)
    return p


print(f"{penality()=}")

model.write("C:/Users/asus/Desktop/ETP_problem/sol/ETP_Problem_test.sol")

with open("C:/Users/asus/Desktop/ETP_problem/sol/ETP_Problem_test.sol", 'w') as f:
    f.write(f"{check_feasible()=}")
    f.write(f"\n{penality()=}")
    f.write(f"\ntimeslot_of_exam={timeslot_of_exam[1:]}")
    f.write("\n")
    f.write(f"exams_at_timeslot={exams_at_timeslot[1:]}")

# Closing the model and environment
model.dispose()
env.dispose()
