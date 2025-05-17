
package lipid;

import adduct.*;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IonizationMode ionizationMode;
    private String adduct;
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;
    private final int PPMTOLERANCE=10;



    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        this.adduct=identifyAdduct();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IonizationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /*public String identifyAdduct() {

        String _adduct = null;

        Map<String, Double> adductMap = ionizationMode == IoniationMode.POSITIVE
                ? AdductList.MAPMZPOSITIVEADDUCTS
                : AdductList.MAPMZNEGATIVEADDUCTS;
                for (String adduct1 : adductMap.keySet()) {
                    for (String adduct2 : adductMap.keySet()) {
                        if (adduct1.equals(adduct2)) continue;

                        for (Peak p1 : groupedSignals) {
                            for (Peak p2 : groupedSignals) {
                                if (p1.equals(p2)) continue;

                                double mz1 = p1.getMz();
                                double mz2 = p2.getMz();

                                Double mass1 = getMonoisotopicMassFromMZ(mz1, adduct1, ionizationMode);
                                Double mass2 = getMonoisotopicMassFromMZ(mz2, adduct2, ionizationMode);


                                if (mass1 != null && mass2 != null) {
                                    int ppmDiffBetweenMasses = Adduct.calculatePPMIncrement(mass1, mass2);
                                    if (ppmDiffBetweenMasses <= PPMTOLERANCE) {
                                        int ppmToMz1 = Adduct.calculatePPMIncrement(mz1, this.mz);
                                        int ppmToMz2 = Adduct.calculatePPMIncrement(mz2, this.mz);

                                        if (ppmToMz1 <= PPMTOLERANCE) {
                                            _adduct = adduct1;
                                            return this.adduct = _adduct;

                                        } else if (ppmToMz2 <= PPMTOLERANCE) {
                                            _adduct = adduct2;
                                            return this.adduct = _adduct;

                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            return this.adduct = _adduct;

    }*/
    public String identifyAdduct() {
        String detectedAdduct = null;

        Map<String, Double> adductMap = ionizationMode == IonizationMode.POSITIVE
                ? AdductList.MAPMZPOSITIVEADDUCTS
                : AdductList.MAPMZNEGATIVEADDUCTS;

        for (String adduct1 : adductMap.keySet()) {
            for (String adduct2 : adductMap.keySet()) {
                if (adduct1.equals(adduct2)) continue;

                for (Peak p1 : groupedSignals) {
                    for (Peak p2 : groupedSignals) {
                        if (p1.equals(p2)) continue;

                        double mz1 = p1.getMz();
                        double mz2 = p2.getMz();

                        Double mass1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1, ionizationMode);
                        Double mass2 = Adduct.getMonoisotopicMassFromMZ(mz2, adduct2, ionizationMode);

                        if (mass1 != null && mass2 != null) {
                            int ppmDiffBetweenMasses = Adduct.calculatePPMIncrement(mass1, mass2);

                            if (ppmDiffBetweenMasses <= PPMTOLERANCE) {
                                int ppmToMz1 = Adduct.calculatePPMIncrement(mz1, this.mz);
                                int ppmToMz2 = Adduct.calculatePPMIncrement(mz2, this.mz);

                                if (ppmToMz1 <= PPMTOLERANCE) {
                                    detectedAdduct = adduct1;
                                    break;
                                } else if (ppmToMz2 <= PPMTOLERANCE) {
                                    detectedAdduct = adduct2;
                                    break;
                                }
                            }
                        }
                    }
                    if (detectedAdduct != null) break;
                }
                if (detectedAdduct != null) break;
            }
            if (detectedAdduct != null) break;
        }

        // 2. Fase 2: si no se detectó nada, probar directamente con this.mz contra cada aducto
        if (detectedAdduct == null) {
            for (String adduct : adductMap.keySet()) {
                Double monoMass = Adduct.getMonoisotopicMassFromMZ(this.mz, adduct, ionizationMode);
                if (monoMass != null) {
                    Double expectedMz = Adduct.getMZFromMonoisotopicMass(monoMass, adduct, ionizationMode);
                    if (expectedMz != null) {
                        int ppm = Adduct.calculatePPMIncrement(this.mz, expectedMz);
                        if (ppm <= PPMTOLERANCE) {
                            detectedAdduct = adduct;
                            break;
                        }
                    }
                }
            }
        }

        this.adduct = detectedAdduct;
        return this.adduct;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {

        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    // !!TODO Detect the adduct with an algorithm or with drools, up to the user.
}
