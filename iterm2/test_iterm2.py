#!/usr/bin/env python3
 # -*- coding: utf-8 -*-

import iterm2

async def main(connection):
    @iterm2.TitleProviderRPC
    async def badge_title(
        badge=iterm2.Reference("badge?"),
        auto_name=iterm2.Reference("autoName?")):
        if badge and auto_name:
            return auto_name + u" \u2014 " + badge
        elif auto_name:
            return auto_name
        elif badge:
            return badge
        else:
            return "Shell"
    await badge_title.async_register(connection, "Name + Badge", "com.iterm2.example.name-and-badge")

iterm2.run_forever(main)
