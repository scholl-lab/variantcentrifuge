package org.apache.commons.lang;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

class ArrayUtilsTest {
    @Test
    void nullInputReturnsNull() {
        assertNull(ArrayUtils.toPrimitive(null));
    }

    @Test
    void emptyArrayReturnsEmptyArray() {
        assertArrayEquals(new int[0], ArrayUtils.toPrimitive(new Integer[0]));
    }

    @Test
    void integerArrayConvertsExactly() {
        assertArrayEquals(new int[] {7, 0, -4}, ArrayUtils.toPrimitive(new Integer[] {7, 0, -4}));
    }

    @Test
    void nullElementRetainsNullPointerException() {
        assertThrows(
                NullPointerException.class,
                () -> ArrayUtils.toPrimitive(new Integer[] {1, null, 3}));
    }
}
