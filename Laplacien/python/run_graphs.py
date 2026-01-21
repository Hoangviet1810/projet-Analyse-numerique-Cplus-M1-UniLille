import subprocess

scripts = ["graph.py", "graph2.py", "graph3.py", "graph4.py"]
processes = []

for s in scripts:
    p = subprocess.Popen(["python3", s])
    processes.append(p)

for p in processes:
    p.wait()
