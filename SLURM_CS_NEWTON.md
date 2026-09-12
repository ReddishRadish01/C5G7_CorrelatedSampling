# GPU SLURM에서 Newton-CS critical rod search 실행

이 프로젝트의 SLURM 기본값은 현재 Windows 로컬 실행 설정과 맞춰져 있다. 서버 GPU에 맞춰 CUDA architecture만 `native`로 자동 선택한다.

## 현재 기본 입력값

| 항목 | 값 |
|---|---:|
| Particles per cycle | 2,000,000 |
| Inactive / active cycles | 100 / 200 |
| Initial rod depth | 64.24 cm |
| Search bounds | 21.42–128.52 cm |
| CS probe / minimum probe | 0.25 / 0.03125 cm |
| CS generations | 10 |
| Newton damping / maximum step | 1.0 / 30.0 cm |
| Maximum Newton iterations | 8 |
| Base seed | 92235922383 |
| Fast math / device LTO | enabled / enabled |

Monte Carlo 및 Newton 입력값은 `scripts/critical_rod_newton_cs.env`에서 바꾼다. GPU, CPU, memory, wall time 등 SLURM 자원은 `scripts/run_critical_rod_newton_cs_slurm.sbatch` 상단의 `#SBATCH` 줄에서 바꾼다.

## 방법 A: Windows에서 자동 업로드 및 제출

프로젝트 루트의 PowerShell에서 다음을 실행한다.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\upload_and_submit_newton_cs.ps1 `
  -Server USER@LOGIN_SERVER `
  -RemoteDir C5G7_runs
```

클러스터에서 CUDA를 module로 불러와야 한다면 `-CudaModule cuda/12.8`처럼 서버의 실제 module 이름을 추가한다. 필요한 경우 `-Partition gpu`, `-Account ACCOUNT_NAME`, `-Port 22`도 지정할 수 있다. 이 명령은 시간표시가 붙은 별도 디렉터리를 만들고, 소스를 전송한 다음 `sbatch`까지 수행한다.

업로드만 하고 제출하지 않으려면 `-UploadOnly`를 추가한다.

## 방법 B: 고정 압축파일을 만들어 수동 전송

먼저 전송용 압축파일을 만든다.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\prepare_slurm_newton_cs_package.ps1
```

생성 파일은 `dist/C5G7_CS_Newton_SLURM_current.tar.gz`이다. 그다음 서버로 전송한다.

```powershell
scp .\dist\C5G7_CS_Newton_SLURM_current.tar.gz USER@LOGIN_SERVER:~/
ssh USER@LOGIN_SERVER
```

서버에서 압축을 풀고 제출한다.

```bash
mkdir -p ~/C5G7_CS_Newton
tar -xzf ~/C5G7_CS_Newton_SLURM_current.tar.gz -C ~/C5G7_CS_Newton
cd ~/C5G7_CS_Newton
mkdir -p results
sbatch scripts/run_critical_rod_newton_cs_slurm.sbatch
```

일회성으로 값만 덮어쓸 수도 있다.

```bash
sbatch scripts/run_critical_rod_newton_cs_slurm.sbatch \
  PARTICLES=4000000 INACTIVE=150 ACTIVE=300 INITIAL_DEPTH_CM=64.24
```

## 컴파일과 실행

기본값 `NO_BUILD=0`이므로 SLURM job 시작 시 `nvcc`로 `build/c5g7`을 컴파일한 뒤 같은 job 안에서 계산한다. 이미 동일한 소스로 컴파일된 실행파일을 재사용할 때만 `NO_BUILD=1`을 준다.

현재 빌드는 `-O3`, `-DNDEBUG`, device LTO, `--use_fast_math`를 사용한다. `CUDA_ARCH=native`이므로 할당된 서버 GPU의 architecture를 사용한다. CUDA module을 쓰는 서버에서는 제출 전에 module을 직접 load하거나 자동 업로더의 `-CudaModule`을 사용한다.

모든 계산 seed는 실행 직전에 48-bit 범위로 맞춘 뒤, 짝수이면 1을 더해서 홀수로 강제한다.

## 진행 확인 및 정지

`JOB_ID`는 `sbatch` 출력의 숫자이다.

```bash
squeue -j JOB_ID
tail -f results/slurm-c5g7-newton-cs-JOB_ID.out
tail -f results/slurm-c5g7-newton-cs-JOB_ID.err
```

계산을 정지하려면 다음을 실행한다.

```bash
scancel JOB_ID
```

세부 iteration 로그, `newton_history.csv`, 최종 `result.txt`는 `results/critical_rod_newton_cs_slurm_JOB_ID_*` 디렉터리에 저장된다.
