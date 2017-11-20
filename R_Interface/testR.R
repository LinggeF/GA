
args=commandArgs(TRUE);

n_modules=args[1]
ga_times=args[2]
ls_times=args[3]
input_path=args[4]
output_path=args[5]

print("n_modules:")
print(n_modules)
print("ga_times:")
print(ga_times)
print("ls_times:")
print(ls_times)
print("input path:")
print(input_path)
print("output_path")
print(output_path)

lines=readLines(input_path)
os=file(output_path,"w")
write(n_modules,os,append = T)
write(ga_times,os,append = T)
write(ls_times,os,append = T)
write(lines,os,append = T)
close(os)
