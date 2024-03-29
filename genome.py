#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#! BioInfo-Genome
#!
#! Copyright (c) 2022, ThomasByr.
#! AGPL-3.0-or-later (https://www.gnu.org/licenses/agpl-3.0.en.html)
#! All rights reserved.
#!
#! Redistribution and use in source and binary forms, with or without
#! modification, are permitted provided that the following conditions are met:
#!
#! * Redistributions of source code must retain the above copyright notice,
#!   this list of conditions and the following disclaimer.
#!
#! * Redistributions in binary form must reproduce the above copyright notice,
#!   this list of conditions and the following disclaimer in the documentation
#!   and/or other materials provided with the distribution.
#!
#! * Neither the name of this software's authors nor the names of its
#!   contributors may be used to endorse or promote products derived from
#!   this software without specific prior written permission.
#!
#! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#! POSSIBILITY OF SUCH DAMAGE.

import sys

if sys.version_info < (3, 10, 6):
    print("Python 3.10.6 or higher is required.")
    sys.exit(1)


if __name__ == "__main__":
    import os
    import dotenv
    import logging
    from src import Tree, App
    from src.helper.logger import init_logger

    dotenv.load_dotenv()
    DEBUG = os.getenv("DEBUG", "false").lower() in {"true", "yes", "1"}
    REBUILD = os.getenv("REBUILD")

    init_logger(logging.DEBUG if DEBUG else logging.INFO)

    overview = Tree(should_rebuild=REBUILD)
    overview.build()

    app = App(overview)
    app.mainloop()
    sys.exit(0)
