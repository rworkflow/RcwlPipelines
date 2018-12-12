cwlVersion: v1.0
class: CommandLineTool
label: Markfiles key-merge reduce
baseCommand: ["python", "/opt/key-merge.py","reduce"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/markfiles-tool:latest

inputs:
  mappedMAF:
    type: File
    inputBinding:
      position: 1
  origMAF:
    type: File
    inputBinding:
      prefix: --maf
  outMAF:
    type: string
    default: reduced.maf
    inputBinding:
      prefix: --output
  appendFlag:
    type: boolean
    default: false
    inputBinding:
      prefix: --append

outputs:
  reducedMAF:
    type: File
    outputBinding:
      glob: $(inputs.outMAF)
