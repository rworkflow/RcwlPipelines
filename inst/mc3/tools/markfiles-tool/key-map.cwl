cwlVersion: v1.0
class: CommandLineTool
label: Markfiles key-merge map
baseCommand: ["python", "/opt/key-merge.py","map"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/markfiles-tool:latest

inputs:
  markList:
    type: File
    inputBinding:
      position: 1
  mergedMAF:
    type: File
    inputBinding:
      prefix: --maf
  useNorm:
    type: int
    default: 0
    inputBinding:
      prefix: --type
  outMAF:
    type: string
    default: mapped.maf
    inputBinding:
      prefix: --output

outputs:
  mappedMAF:
    type: File
    outputBinding:
      glob: $(inputs.outMAF)
