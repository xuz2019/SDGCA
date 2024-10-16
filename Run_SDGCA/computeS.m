function ML = computeS(NWCA, MLA)
    ML = NWCA;
    ML(MLA == 0) = 0;
    ML = ML - min(min(ML));
    ML = ML / max(max(ML));
    ML = ML / 5 + 0.8;
    ML(ML == 0.8) = 0;
end

