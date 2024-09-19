import numpy as np

from primalbedtools.bedfiles import BedLine


def create_mapping_array(
    msa: dict[str, str], from_id: str, to_id: str
) -> tuple[np.ndarray, dict[int, int]]:
    """
    Returns a mapping array and a dict of from_genome index to msa index
    from_id will be array[0, :] and id_to will be array[1, :]
    """
    # Check if IDs are in the MSA
    if from_id not in msa:
        raise ValueError(f"ID {from_id} not found in ({', '.join(msa.keys())})")
    if to_id not in msa:
        raise ValueError(f"ID {to_id} not found in ({', '.join(msa.keys())})")

    # Check for same names
    if from_id == to_id:
        raise ValueError("IDs are the same")

    # Check for different lengths
    if len(msa[from_id]) != len(msa[to_id]):
        raise ValueError("MSA lengths are different")

    # The +1 is needed to account for edge case with primer at the end, due to non-inclusive slicing
    msa_to_genome = np.full([2, len(msa[from_id]) + 1], None)  # type: ignore

    # populate with from_genome indexes
    from_seq = msa[from_id]
    from_index = 0
    for msa_index in range(len(from_seq)):
        if from_seq[msa_index] not in {"", "-"}:
            msa_to_genome[0, msa_index] = from_index
            from_index += 1
    msa_to_genome[0, -1] = from_index

    # to genome indexes
    to_seq = msa[to_id]
    to_index = 0
    for msa_index in range(len(to_seq)):
        if to_seq[msa_index] not in {"", "-"}:
            msa_to_genome[1, msa_index] = to_index
            to_index += 1
    msa_to_genome[1, -1] = to_index

    # Create a dict of primary ref to msa
    from_index_to_msa_index = {}
    for msa_index, from_index in enumerate(msa_to_genome[0]):
        if from_index is not None:
            from_index_to_msa_index[from_index] = msa_index

    return msa_to_genome, from_index_to_msa_index


def remap(
    from_id: str,
    to_id: str,
    bedlines: list[BedLine],
    msa: dict[str, str],
):
    msa_to_genome, from_index_to_msa_index = create_mapping_array(msa, from_id, to_id)

    for bedline in bedlines:
        # Guard for bedlines to other chromosomes
        if bedline.chrom != from_id:
            continue

        msa_start = from_index_to_msa_index[bedline.start]
        msa_end = from_index_to_msa_index[bedline.end]

        # Check for perfect mapping
        if None not in msa_to_genome[:, msa_start:msa_end]:
            bedline.start = msa_to_genome[1, msa_start]
            bedline.end = msa_to_genome[1, msa_end - 1] + 1
            bedline.chrom = to_id
            continue

        # Check for primer not in the new reference
        if np.flatnonzero(msa_to_genome[1, msa_start:msa_end]).size == 0:
            print(f"{bedline.primername} not found in new reference")
            # revert to original
            continue

        # Handle non 3' gaps
        new_ref_slice = msa_to_genome[1, msa_start:msa_end]

        if (
            new_ref_slice[-1] if bedline.strand == "+" else new_ref_slice[0]
        ) is not None:
            if bedline.strand == "+":
                bedline.end = msa_to_genome[1, msa_end - 1] + 1
                bedline.start = max(bedline.end - len(bedline.sequence), 0)
            else:
                bedline.start = msa_to_genome[1, msa_start]
                bedline.end = min(
                    bedline.start + len(bedline.sequence), len(msa_to_genome[1])
                )
            bedline.chrom = to_id
            continue
        else:
            print(f"{bedline.primername} 3' gap found in new reference")

        # Handle 3' gaps
        # At this point at least one base is 'mapped'
        # Find the next valid 3' base
        if bedline.strand == "+":
            for i in range(msa_end, msa_to_genome.shape[1]):
                if msa_to_genome[1, i] is not None:
                    bedline.end = msa_to_genome[1, i] + 1
                    bedline.start = max(bedline.end - len(bedline.sequence), 0)
                    bedline.chrom = to_id
                    break
            continue
        else:
            for i in range(msa_start, -1, -1):
                print(i)
                if msa_to_genome[1, i] is not None:
                    bedline.start = msa_to_genome[1, i]
                    bedline.end = min(
                        bedline.start + len(bedline.sequence), len(msa_to_genome[1])
                    )
                    bedline.chrom = to_id
                    break
            continue

    return bedlines
