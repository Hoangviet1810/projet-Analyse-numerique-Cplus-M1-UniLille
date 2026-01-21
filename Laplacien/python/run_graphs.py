import subprocess
import os
import glob

# thư mục chứa run_graphs.py và các script .py
python_dir = os.path.dirname(os.path.abspath(__file__))

# danh sách script
scripts = ["graph.py", "graph2.py", "graph3.py", "graph4.py"]
processes = []

# chạy tất cả script song song
for s in scripts:
    path = os.path.join(python_dir, s)
    p = subprocess.Popen(["python3", path])
    processes.append(p)

# đợi tất cả kết thúc
for p in processes:
    p.wait()

# --- xoá các file kết quả .txt trong thư mục python ---
txt_files = glob.glob(os.path.join(python_dir, "residus_*.txt"))
txt_files.append(os.path.join(python_dir, "T_final.txt"))

for f in txt_files:
    if os.path.exists(f):
        os.remove(f)
        print(f"Supprimé: {f}")

print("Tous les fichiers de résultats ont été supprimés.")