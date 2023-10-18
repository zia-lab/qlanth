
BeginPackage["Slurmy`"];

(*
For this to work the resources need to have been already allocated in the cluster.
For this purpose the command below may be used.
I probably can request more resources but presently this seems like more than plenty.
*)

(*
mem=32
numCpus=10
hours=24
for machname in alpha beta gamma epsilon iota upsilon
do
echo "$machname"
salloc -J $machname -N 1-1 -n $numCpus --time=${hours}:00:00 --mem=${mem}g -p batch srun --pty bash &
done
*)

goodnodes = {"alpha", "beta", "gamma", "epsilon", "iota", "upsilon"};

SlurmQueue::usage = "
SlurmQueue[] returns a Dataset representing the current jobs queued 
in the SLURM workload manager for the current user. 

The dataset includes columns:
  - 'JOBID': the unique job identifier.
  - 'NAME': the name of the job.
  - 'STATE': the current state of the job (e.g., RUNNING, PENDING).
  - 'NODES': the number of nodes used by the job.
  - 'NODELIST(REASON)': the names of the nodes where the job is running.
  - 'CPUS': the number of CPUs allocated for the job.
  - 'MIN_MEMORY': the amount of memory requested by the job in gigabytes.

The function queries the SLURM workload manager using the 'squeue' command 
and requires that 'squeue' be available on the system path.

Example:
  ds = SlurmQueue[];
  ds[All, 'JOBID'] (* Lists all job IDs *)

Notes:
  - The function assumes that the SLURM environment is properly set up and that 
    the user has permission to query job details.
  - The exact details and formatting of the returned dataset might vary based 
    on the SLURM configuration and version.
";

SlurmQueue[] := 
 Module[{cmd, result, lines, headers, data}, 
  cmd = {"squeue", "-u", $UserName, 
    "--format=%i %j %T %M %D %R %C %m"};
  result = RunProcess[cmd];
  If[result["ExitCode"] != 0, Return[$Failed, Module]];
  lines = StringSplit[result["StandardOutput"], "\n"];
  headers = StringTrim /@ StringSplit[First@lines];
  data = Rest[lines];
  (*Parsing the rest of the data*)
  data = Map[
    StringTrim /@ 
      StringSplit[#, RegularExpression["(?<=\\S)\\s+(?=\\S)"]] &, 
    data];
  dataslurm = Dataset[AssociationThread[headers -> #] & /@ data];
  Return[dataslurm]]

SlurmLauncher::usage = "
SlurmLauncher[] launches Mathematica kernels on the compute nodes where the current user 
has jobs running under the SLURM workload manager, provided these nodes are members of a predefined 'goodnodes' list.

The function does the following:
  - Queries the current SLURM queue for jobs using the SlurmQueue[] function.
  - Filters jobs that are running on nodes listed in the 'goodnodes' list.
  - For each of these filtered jobs, it launches a Mathematica kernel on the corresponding node using the number of allocated CPUs.
  - Uses a specific KernelCommand path to launch the Mathematica kernel.

Returns:
  - A list containing three elements:
    1. 'totalKernels': The total number of Mathematica kernels launched.
    2. 'totalMemory': The total amount of memory (in gigabytes) associated with the launched kernels.
    3. 'runningNodes': A list of nodes, CPUs, and memory where the kernels were launched. Each element of this list is of the form {nodeName, numCPUs, memory}.

Assumptions and Requirements:
  - The function assumes that the SLURM environment is properly set up and that the user has permission to query job details.
  - The 'goodnodes' list must be predefined and contains the names of nodes considered suitable for kernel launch.
  - The function uses a specific path in the 'KernelCommand' option to launch Mathematica, which should be adjusted based on the system configuration.
  - Requires 'squeue' to be available on the system path.

Example:
  result = SlurmLauncher[];
  Print[result[[1]]]; (* Prints the total number of kernels launched *)
";

SlurmLauncher[] := (
  slurm = SlurmQueue[];
  runningNodes = 
   Select[Normal[slurm], MemberQ[goodnodes, #[["NAME"]]] &];
  runningNodes = {#[["NODELIST(REASON)"]], 
      ToExpression[#[["CPUS"]]], ToExpression[StringReplace[#[["MIN_MEMORY"]],"G"->""]]} & /@ runningNodes;
  totalKernels = 0;
  totalMemory = 0;
  Do[(
    {host, numCPUs, mem} = runny;
    kernelSSH = "ssh://" <> host;
    Print["Launching ", numCPUs, " kernels on ", host];
    totalKernels += numCPUs;
    totalMemory += mem;
    LaunchKernels[KernelConfiguration[kernelSSH,
      "KernelCommand" -> "/gpfs/runtime/opt/mathematica/13.2.0/bin/MathKernel", 
      "KernelCount" -> numCPUs
      ]]
    ),
    {runny, runningNodes}];
  Return[{totalKernels,totalMemory,runningNodes}]
  )

EndPackage[];

If[
  SameQ[$FrontEnd, Null],
  (
  {numKernels, totalMemory, theNodes} = SlurmLauncher[];
  Print[{numKernels, totalMemory}];
  ParallelTable[
      (identity = {$MachineName, $ProcessID};
      Print["Executing on ", identity, " for i = ", i];
      ),
      {i, 1, numKernels}
  ];
  )
]
