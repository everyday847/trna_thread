import typing

class Sequence:
    def __init__(self, sequence: str, secstruct=
        ".(((((((..((((...........)))).(((((.......)))))........................(((((.......))))))))))))....") -> None:
        self.sequence: str = sequence
        self.secstruct: str = secstruct

    def __repr__(self) -> str:
        return "SEQUENCE: {}\nSECSTRUCT: {}\n".format(self.sequence, self.secstruct)

    def __len__(self) -> int: return len(self.sequence)
