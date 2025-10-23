import dataclasses

UPREGULATED = "upregulated"
DOWNREGULATED = "downregulated"
NOT_REGULATED = "not_regulated"

@dataclasses.dataclass
class GeneNode:
    name: str
    fold_change: float
    p_value: float
    padj: float
    combined_score: float

    def is_regulated(self) -> str:
        if self.fold_change > 0:
            return UPREGULATED
        elif self.fold_change < 0:
            return DOWNREGULATED
        else:
            return NOT_REGULATED
        
    
    @classmethod
    def from_dict(cls, data: dict) -> "GeneNode":
        """Create a GeneNode from a dictionary (handles naming differences too)."""
        return cls(
            name = data["name"],
            fold_change = data["log2FoldChange"] if "log2FoldChange" in data else 0.0,
            p_value = data["pvalue"] if "pvalue" in data else 0.0,
            padj = data["padj"] if "padj" in data else 0.0,
            combined_score = data["combined_score"] if "combined_score" in data else 0.0,
        )
    
    def to_dict(self):
        """Serialize GeneNode to a dictionary with only basic types"""
        return {
            "name": self.name,
            "fold_change": self.fold_change,
            "p_value": self.p_value,
            "padj": self.padj,
            "combined_score": self.combined_score,
        }