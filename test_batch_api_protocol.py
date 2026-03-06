from polymerase_tm.batch import from_csv, to_csv

print("Testing PROTOCOL batch...")
results = from_csv("test_input_protocol.csv", action="protocol", name_col="name", fwd_col="fwd", rev_col="rev", template_col="template")
to_csv(results, "test_output_protocol.csv")
print("PROTOCOL done.")

with open("test_output_protocol.csv", encoding="utf-8") as f:
    print(f.read()[:500])
