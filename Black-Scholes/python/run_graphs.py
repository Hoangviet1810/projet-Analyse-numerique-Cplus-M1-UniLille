import subprocess
import os
import glob

# thư mục chứa run_graphs.py và các script .py
python_dir = os.path.dirname(os.path.abspath(__file__))

# danh sách script
scripts = ["graph1.py", "graph2.py"]
processes = []

# chạy tất cả script song song
for s in scripts:
    path = os.path.join(python_dir, s)
    p = subprocess.Popen(["python3", path])
    processes.append(p)

# đợi tất cả kết thúc
for p in processes:
    p.wait()

# --- xoà các file kết quả .txt trong thư mục python ---
txt_files = glob.glob(os.path.join(python_dir, "BS_implicit.txt"))
txt_files.append(os.path.join(python_dir, "BS_implicit.txt"))

for f in txt_files:
    if os.path.exists(f):
        os.remove(f)
        print(f"Supprimé: {f}")

print("Tous les fichiers de résultats ont été supprimés.")