from polymerase_tm.batch import from_csv, to_csv

print("Testing SDM batch...")
results = from_csv("test_input_sdm.csv", action="sdm", template_col="template", mutation_col="mutation", name_col="name")
to_csv(results, "test_output_sdm.csv")
print("SDM done.")

with open("test_output_sdm.csv") as f:
    print(f.read()[:500])
