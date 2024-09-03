package test.util;

import main.java.util.MGFUtils;
import org.junit.jupiter.api.Test;
import java.util.HashMap;

public class MGFUtilsTest {

    @Test
    void testGet_mod_from_peptide() {
        String seq = "-17.027SGLQENAFVNMKPSQILQ+0.984TVK";
        HashMap<Integer,Double> map = MGFUtils.get_mod_from_peptide(seq);
        org.junit.jupiter.api.Assertions.assertEquals(map.get(0), -17.027);
        org.junit.jupiter.api.Assertions.assertEquals(map.get(18), 0.984);
    }
}
