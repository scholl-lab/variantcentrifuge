package org.apache.commons.lang;

/** Minimal binary compatibility for mmtf-codec 1.0.10. */
public final class ArrayUtils {
    private ArrayUtils() {}

    public static int[] toPrimitive(Integer[] array) {
        if (array == null) {
            return null;
        }
        int[] result = new int[array.length];
        for (int index = 0; index < array.length; index++) {
            result[index] = array[index].intValue();
        }
        return result;
    }
}
