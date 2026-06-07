import json
import traceback
import sys

def verify():
    notebook_path = "pca_analysis.ipynb"
    print(f"Reading notebook: {notebook_path}")
    
    with open(notebook_path, "r", encoding="utf-8") as f:
        nb = json.load(f)
    
    # Extract code cells
    code_cells = []
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            source = "".join(cell.get("source", []))
            code_cells.append(source)
            
    print(f"Found {len(code_cells)} code cells.")
    
    # Combine code cells, stripping IPython magic lines (like %matplotlib)
    clean_code_lines = []
    for i, cell_code in enumerate(code_cells):
        clean_code_lines.append(f"# --- Cell {i+1} ---")
        for line in cell_code.splitlines():
            if line.strip().startswith("%"):
                # Comment out magic commands
                clean_code_lines.append(f"# {line}")
            else:
                clean_code_lines.append(line)
                
    full_code = "\n".join(clean_code_lines)
    
    # Execute the code in a clean local namespace
    locs = {}
    print("Executing notebook code...")
    try:
        exec(full_code, globals(), locs)
        print("\nSUCCESS: All cells executed successfully without errors!")
        print("SUCCESS: All biological invariants (assertions) passed!")
        return True
    except Exception as e:
        print("\nERROR occurred during notebook execution:")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = verify()
    sys.exit(0 if success else 1)
