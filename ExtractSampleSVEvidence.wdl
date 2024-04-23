version 1.0

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

# Extract sample-level evidence from merged batch evidence
workflow ExtractSampleSVEvidence {
  input {
    String batch_id
    String batch_bucket_path
    Array[String] sample_ids
    String gatk_docker
    File rd_header
    RuntimeAttr? runtime_attr_override
  }

  String batch_bucket_path_clean = sub(batch_bucket_path, "/*$", "")
  File batch_pe_file = batch_bucket_path_clean + "/" + batch_id + ".PE.txt.gz"
  File batch_sr_file = batch_bucket_path_clean + "/" + batch_id + ".SR.txt.gz"
  File batch_rd_file = batch_bucket_path_clean + "/" + batch_id + ".RD.txt.gz"
  File batch_baf_file = batch_bucket_path_clean + "/" + batch_id + ".BAF.txt.gz"

  call GetPE {
    input:
      batch_pe_file = batch_pe_file,
      sample_ids = sample_ids,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  call GetSR {
    input:
      batch_sr_file = batch_sr_file,
      sample_ids = sample_ids,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  call GetRD {
    input:
      batch_rd_file = batch_rd_file,
      rd_header = rd_header,
      sample_ids = sample_ids,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  call GetBAF {
    input:
      batch_baf_file = batch_baf_file,
      sample_ids = sample_ids,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    Array[File] pe_files = GetPE.pe_files
    Array[File] sr_files = GetSR.sr_files
    Array[File] rd_files = GetRD.rd_files
    Array[File] baf_files = GetBAF.baf_files
  }
}

task GetPE {
  input {
    File batch_pe_file
    Array[String] sample_ids
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float batch_pe_file_size = size(batch_pe_file, "GiB")
  Int vm_disk_size = ceil(batch_pe_file_size * 3)
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])  
  File sample_ids_file = write_lines(sample_ids)

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir 'PE'
    gzip --decompress --stdout '~{batch_pe_file}' \
      | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($7 in a){print > ("PE/" $7 ".PE.txt")}' '~{sample_ids_file}' -
    find 'PE' -name '*.PE.txt' -type f -exec bgzip '{}' \; 
  >>>

  output {
    Array[File] pe_files = glob("PE/*.PE.txt.gz")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetSR {
  input {
    File batch_sr_file
    Array[String] sample_ids
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float batch_sr_file_size = size(batch_sr_file, "GiB")
  Int vm_disk_size = ceil(batch_sr_file_size * 3)
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])  
  File sample_ids_file = write_lines(sample_ids)

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir 'SR'
    gzip --decompress --stdout '~{batch_sr_file}' \
      | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($5 in a){print > ("SR/" $5 ".SR.txt")}' '~{sample_ids_file}' -
    find 'SR' -name '*.SR.txt' -type f -exec bgzip '{}' \; 
  >>>

  output {
    Array[File] sr_files = glob("SR/*.SR.txt.gz")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetRD {
  input {
    File batch_rd_file
    File rd_header
    Array[String] sample_ids
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float batch_rd_file_size = size(batch_rd_file, "GiB")
  Int vm_disk_size = ceil(batch_rd_file_size * 5)
  RuntimeAttr default_attr  = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])  
  File sample_ids_file = write_lines(sample_ids)

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir 'RD'
    gzip --decompress --stdout '~{batch_rd_file}' \
    | awk -F'\t' 'NR==FNR{a[$1]; next}
        FNR==1{for(i=1;i<=NF;++i)b[$i]=i}
        FNR>1{for(id in b){if(id in a) print $1,$2,$3,$(b[id]) > ("RD/" id ".RD.txt")}}' OFS="\t" '~{sample_ids_file}' -
    while IFS= read -r file; do
      bn="$(basename "${file}")"
      sample_id="${bn%.RD.txt}"
      printf '@RG\tID:GATKCopyNumber\tSM:%s\n' "${sample_id}" \
        | cat '~{rd_header}' - "${file}" > "RD/${sample_id}.tmp"
      mv "RD/${sample_id}.tmp" "${file}"
      bgzip "${file}"
    done < <(find 'RD' -name '*.RD.txt' -type f -print)
  >>>

  output {
    Array[File] rd_files = glob("RD/*.RD.txt.gz")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetBAF {
  input {
    File batch_baf_file
    Array[String] sample_ids
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float batch_baf_file_size = size(batch_baf_file, "GiB")
  Int vm_disk_size = ceil(batch_baf_file_size * 3)
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])  
  File sample_ids_file = write_lines(sample_ids)

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir 'BAF'
    gzip --decompress --stdout '~{batch_baf_file}' \
      | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($4 in a){print > ("BAF/" $4 ".BAF.txt")}' '~{sample_ids_file}' -
    find 'BAF' -name '*.BAF.txt' -type f -exec bgzip '{}' \; 
  >>>

  output {
    Array[File] baf_files = glob("BAF/*.BAF.txt.gz")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
