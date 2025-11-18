from __future__ import annotations

from typing import Dict, List


COLOR_PRESETS: Dict[str, List[str]] = {
    "Basic": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"],
    "Steve": ["#F7FFA1", "#756027", "#45FF42", "#FF9100", "#1244D4", "#FF0000", "#1F1F1F", "#FFA1FC", "#9903A1", "#BBFAF5", "#0E661C"],
    "Pastel": ["#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"],
    "Rainbow": ["#E6194B", "#3CB44B", "#FFE119", "#0082C8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", "#D2F53C", "#FABEBE"],
    "Warm": ["#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"],
    "Earth": ["#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"],
    "Vibrant": ["#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080"],
    "Light": ["#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD"],
    "Neon": ["#E6FB04", "#00FF00", "#FF00CC", "#9D00FF", "#00FFFF", "#FF0000"],
    "Dark": ["#8A2BE2", "#000000", "#008B8B", "#006400", "#483D8B", "#B8C709"],
    "Jade": ["#94E8B4", "#72BDA3", "#5E8C61", "#4E6151", "#3B322C", "#CCAD99"],
    "Beach ball": ["#25CED1", "#FCEADE", "#FF8A5B", "#EA526F", "#FFFFFF", "#FFC300"],
}


def preset_colors(name: str, n: int) -> List[str]:
    palette = COLOR_PRESETS.get(name) or COLOR_PRESETS["Basic"]
    if len(palette) >= n:
        return palette[:n]
    cycles = (n // len(palette)) + 1
    return (palette * cycles)[:n]
