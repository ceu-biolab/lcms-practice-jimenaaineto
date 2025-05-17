package adduct;

import lipid.IonizationMode;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Adduct {

    public static int extractMultimer(String adduct) {

        Pattern pattern = Pattern.compile("\\[(\\d*)M");
        Matcher matcher = pattern.matcher(adduct);
        if (matcher.find() && !matcher.group(1).isEmpty()) {
            return Integer.parseInt(matcher.group(1));
        }
        return 1;

    }

    public static int extractCharge(String adduct) {
            Pattern pattern = Pattern.compile("(\\d*)([\\+\\-−])\\]");
            Matcher matcher = pattern.matcher(adduct);

       if (matcher.find()) {
            String number = matcher.group(1);
            String sign = matcher.group(2);
            int charge = number.isEmpty() ? 1 : Integer.parseInt(number);
            if (sign.equals("-") || sign.equals("−")) charge *= -1;
            return charge;
        }
        return 1;

    }

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct, IonizationMode ioniation) {
        Double massToSearch;
        if(mz!= null||adduct!=null) {
            int multimer = extractMultimer(adduct);
            int charge = extractCharge(adduct);
            Map<String, Double> adductMap = ioniation == IonizationMode.POSITIVE ? AdductList.MAPMZPOSITIVEADDUCTS : AdductList.MAPMZNEGATIVEADDUCTS;
            massToSearch = adductMap.get(adduct);
            return (mz * charge + massToSearch) / multimer;
        }

            return null;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct, IonizationMode ioniation) {
        Double massToSearch;
        if(monoisotopicMass!= null||adduct!=null) {
            int multimer = extractMultimer(adduct);
            int charge = extractCharge(adduct);
            Map<String, Double> adductMap = ioniation == IonizationMode.POSITIVE ? AdductList.MAPMZPOSITIVEADDUCTS : AdductList.MAPMZNEGATIVEADDUCTS;
            massToSearch = adductMap.get(adduct);
            return (monoisotopicMass * multimer - massToSearch) / charge;
        }

        return null;

    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }




}
