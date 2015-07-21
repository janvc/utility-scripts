# Copyright 2015 Jan von Cosel
#
# This file is part of gauss_viban.
#
# gauss_viban if free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gauss_viban is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have recieved a copy of the GNU General Public License
# along with astronx. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#

FC = gfortran
MKDIR_P = mkdir -p

SRC_DIR = source
OBJ_DIR = objects
MOD_DIR = modules
BIN_DIR = bin

GAUSS_VIBAN = $(addprefix $(BIN_DIR)/,gauss_viban)

CFLAGS = -Wall -Wextra -g --pedantic-errors -std=f2008 -O0 -fall-intrinsics
LFLAGS = -llapack
VPATH = $(SRC_DIR):   \
        $(OBJ_DIR)

OBJECTS = $(addprefix $(OBJ_DIR)/,\
                 routines.o       \
                 gauss_viban.o    \
                 )

all: directories $(GAUSS_VIBAN)

directories: $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

$(GAUSS_VIBAN): $(OBJECTS)
	$(FC) -o $@ $^ $(LFLAGS)

$(OBJ_DIR)/%.o: %.f90
	$(FC) -c $(CFLAGS) -J$(MOD_DIR) $< -o $@

$(OBJ_DIR):
	$(MKDIR_P) $(OBJ_DIR)
  
$(MOD_DIR):
	$(MKDIR_P) $(MOD_DIR)
  
$(BIN_DIR):
	$(MKDIR_P) $(BIN_DIR)

clean:
	test -d $(BIN_DIR) && rm -r $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)
