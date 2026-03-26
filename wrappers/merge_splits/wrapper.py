from pathlib import Path

import h5py


def main():
    output_path = Path(str(snakemake.output[0])).resolve()  # noqa: F821
    input_paths = [Path(path).resolve() for path in snakemake.input]  # noqa: F821

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        output_path.unlink()

    with h5py.File(output_path, "w") as out:
        out.attrs["merge_type"] = "external_links"
        out.attrs["inputs"] = [str(path) for path in input_paths]

        for index, path in enumerate(input_paths, start=1):
            roi_name = f"roi{index:02d}"
            if not path.exists():
                raise FileNotFoundError(str(path))
            out[roi_name] = h5py.ExternalLink(str(path), "/")


main()
