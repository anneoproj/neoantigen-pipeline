import argparse
import subprocess
from pathlib import Path


def run_mixmhc2(input_file, output_file, command, sample_id=None):
    if not command:
        raise ValueError(
            "mixmhc2_command is required when running mixmhc2_only workflow."
        )

    try:
        command = command.format(
            input_file=input_file,
            output_file=output_file,
            sample_id=sample_id or "",
            pairs_file=input_file,
        )
    except KeyError as error:
        raise ValueError(f"Invalid mixmhc2 command template: missing {error}")

    result = subprocess.run(
        command,
        shell=True,
        check=True,
        executable="/bin/bash",
    )
    if result.returncode != 0:
        raise RuntimeError(f"MixMHC2 command failed with exit code {result.returncode}")

    if output_file and not Path(output_file).exists():
        raise FileNotFoundError(f"MixMHC2 output not found: {output_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--sample_id")

    args = parser.parse_args()
    run_mixmhc2(
        input_file=args.input_file,
        output_file=args.output_file,
        command=args.command,
        sample_id=args.sample_id,
    )


if __name__ == "__main__":
    main()
